package ahrd.test;

import static ahrd.controller.Settings.getSettings;
import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.assertNotNull;
import static junit.framework.Assert.assertTrue;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.Vector;

import org.junit.Before;
import org.junit.Test;

import ahrd.controller.Utils;
import ahrd.exception.MissingInterproResultException;
import ahrd.model.BlastResult;
import ahrd.model.DomainScoreCalculator;
import ahrd.model.InterproResult;
import ahrd.model.Protein;

public class DomainScoreCalculatorTest {

	protected static InterproResult interpro1 = new InterproResult("IPR000001",
			"domain one", "Domain", 0.1);
	protected static InterproResult interpro2 = new InterproResult("IPR000002",
			"domain two", "Domain", 0.2);
	protected static InterproResult interpro3 = new InterproResult("IPR000003",
			"domain three", "Domain", 0.3);
	protected static InterproResult interpro4 = new InterproResult("IPR000004",
			"domain four", "Domain", 0.4);
	protected static InterproResult interpro5 = new InterproResult("IPR000005",
			"domain five", "Domain", 0.5);
	protected static InterproResult interpro6 = new InterproResult("IPR000006",
			"domain six", "Domain", 0.6);

	protected void mockBlastResultAccessionsToInterproIds() {
		Map<String, Set<String>> brAccs2iprIds = new HashMap<String, Set<String>>();
		brAccs2iprIds.put(
				"accession_1",
				new HashSet<String>(Arrays.asList(new String[] {
						interpro2.getId(), interpro1.getId() })));
		brAccs2iprIds.put(
				"accession_2",
				new HashSet<String>(Arrays.asList(new String[] {
						interpro1.getId(), interpro3.getId() })));
		brAccs2iprIds.put(
				"accession_3",
				new HashSet<String>(Arrays.asList(new String[] { interpro4
						.getId() })));
		brAccs2iprIds.put(
				"accession_4",
				new HashSet<String>(Arrays.asList(new String[] { interpro5
						.getId() })));
		// accession_5 has no InterproIDs assigned!
		brAccs2iprIds.put(
				"accession_6",
				new HashSet<String>(
						Arrays.asList(new String[] { interpro1.getId(),
								interpro2.getId(), interpro6.getId() })));
		DomainScoreCalculator
				.setBlastResultAccessionsToInterproIds(brAccs2iprIds);
	}

	protected Protein mockProteinWithBlastAndInterpoResults() {
		Protein prot = TestUtils.mockProtein();
		// Interpro-Domains, it has been annotated with:
		prot.getInterproResults().add(interpro1);
		prot.getInterproResults().add(interpro2);
		prot.getInterproResults().add(interpro3);
		// BlastResult found in Uniprot/Swissprot for above Protein,
		// with accessions named 'accession_X' X = 1,2,...,5
		prot.getBlastResults().put("swissprot", TestUtils.mockBlastResults());
		// A single hit in TAIR:
		List<BlastResult> tairHits = new ArrayList<BlastResult>();
		tairHits.add(TestUtils.mockBlastResult("accession_6", 0.001,
				"description six", 1, 20, 1, 20, 200, 100.0, "tair",
				new HashSet<String>()));
		prot.getBlastResults().put("tair", tairHits);
		return prot;
	}

	protected void mockInterproDatabase() {
		Map<String, InterproResult> interproDb = new HashMap<String, InterproResult>();
		interproDb.put(interpro1.getId(), interpro1);
		interproDb.put(interpro2.getId(), interpro2);
		interproDb.put(interpro3.getId(), interpro3);
		interproDb.put(interpro4.getId(), interpro4);
		interproDb.put(interpro5.getId(), interpro5);
		interproDb.put(interpro6.getId(), interpro6);
		InterproResult.setInterproDb(interproDb);
	}

	@Before
	public void initialiseInterproDb() throws IOException {
		TestUtils.initTestSettings();
	}

	@Test
	public void testConstructVectorSpaceModel() {
		// Protein to be tested:
		Protein prot = mockProteinWithBlastAndInterpoResults();
		// Setup memory database 'BlastResultAccessionsToInterproIds'
		mockBlastResultAccessionsToInterproIds();
		// Test the construction of the vector space model:
		SortedSet<String> vectorSpaceModelToTest = DomainScoreCalculator
				.constructVectorSpaceModel(prot);
		SortedSet<String> vectorSpaceModelExpected = new TreeSet<String>(
				Arrays.asList(new String[] { interpro1.getId(),
						interpro2.getId(), interpro3.getId(),
						interpro4.getId(), interpro5.getId(), interpro6.getId() }));
		assertEquals(
				"Construction of vector space model should return the following dimension in their alphabetical order: IPR00000X, X=1,2,...,6",
				vectorSpaceModelExpected, vectorSpaceModelToTest);

	}

	@Test
	public void testConstructDomainWeightVectors()
			throws MissingInterproResultException {
		// Protein to be tested:
		Protein prot = mockProteinWithBlastAndInterpoResults();
		// Setup memory database 'BlastResultAccessionsToInterproIds'
		mockBlastResultAccessionsToInterproIds();
		// Setup memory database 'InterproDb':
		mockInterproDatabase();
		// Construct the domain-weight vectors for the protein and its
		// BlastResults:
		DomainScoreCalculator.constructDomainWeightVectors(prot,
				DomainScoreCalculator.constructVectorSpaceModel(prot));
		// Assure that above vectors have been constructed correctly:
		assertNotNull(prot.getDomainWeights());
		assertEquals(
				"The proteins domain-weight vector is not as expected.",
				new Vector<Double>(Arrays.asList(new Double[] { 0.1, 0.2, 0.3,
						0.0, 0.0, 0.0 })), prot.getDomainWeights());
		// Just test one Swissprot BlastResult, 'accession_1':
		assertNotNull(prot.getBlastResults().get("swissprot").get(0)
				.getDomainWeights());
		assertEquals(
				"The domain-weight vector of BlastResult 'accession_1' is not as expected.",
				new Vector<Double>(Arrays.asList(new Double[] { 0.1, 0.2, 0.0,
						0.0, 0.0, 0.0 })),
				prot.getBlastResults().get("swissprot").get(0)
						.getDomainWeights());
		// ... and the unique TAIR BlastResult, 'accession_6':
		assertNotNull(prot.getBlastResults().get("tair").get(0)
				.getDomainWeights());
		assertEquals(
				"The domain-weight vector of BlastResult 'accession_6' is not as expected.",
				new Vector<Double>(Arrays.asList(new Double[] { 0.1, 0.2, 0.0,
						0.0, 0.0, 0.6 })), prot.getBlastResults().get("tair")
						.get(0).getDomainWeights());
		// BlastResult 'accession_5' should have the ZERO-Vector:
		assertNotNull(prot.getBlastResults().get("swissprot").get(4)
				.getDomainWeights());
		assertEquals(
				"The domain-weight vector of BlastResult 'accession_5' is not as expected.",
				new Vector<Double>(Arrays.asList(new Double[] { 0.0, 0.0, 0.0,
						0.0, 0.0, 0.0 })),
				prot.getBlastResults().get("swissprot").get(4)
						.getDomainWeights());
	}

	@Test
	public void testDomainWeightSimilarity() {
		List<Double> x = Arrays.asList(new Double[] { 0.1, 0.2, 0.3, 0.4, 0.5,
				0.6 });
		List<Double> y = Arrays.asList(new Double[] { 0.2, 0.3, 0.4, 0.5, 0.6,
				0.7 });
		// sim(x,y) = dot-product(x,y) / (||x||*||y||)
		// "In R:" ;-)
		// x = seq(0.1, 0.6, by=0.1)
		// y = seq(0.2, 0.7, by=0.1)
		// sum(x*y)/(sqrt(sum(x^2)) * sqrt(sum(y^2)))
		// = 0.9958408
		Double dws = DomainScoreCalculator.domainWeightSimilarity(x, y);
		assertNotNull(dws);
		assertEquals(0.9958408, dws, 0.000001);

		// Verify that any protein without domain annotation results in zero
		// similarity:
		x = Utils.zeroList(6);
		dws = DomainScoreCalculator.domainWeightSimilarity(x, y);
		assertNotNull(dws);
		assertEquals(0.0, dws, 0.0);
	}

	@Test
	public void testParseBlastResultsDomainAnnotations() throws IOException {
		// Setup
		DomainScoreCalculator
				.setBlastResultAccessionsToInterproIds(new HashMap<String, Set<String>>());
		DomainScoreCalculator
				.setBlastResultAccessionsToPfamIds(new HashMap<String, Set<String>>());
		getSettings().setComputeDomainSimilarityOn("interpro");
		getSettings().setPathToBlastResultsDomainAnnotation(
				"./test/resources/blast_results_domain_annotations_1.tbl");
		// Test
		DomainScoreCalculator.parseBlastResultsDomainAnnotations();
		assertTrue(
				"Should have loaded some Blast Hit Accessions to Domain Annotation Mappings.",
				!DomainScoreCalculator.getBlastResultAccessionsToInterproIds()
						.isEmpty());
	}

	@Test
	public void testParseDomainAnnotationLine() {
		// Setup
		getSettings().setComputeDomainSimilarityOn("interpro");
		String iprAnno1 = "\"SHEEP\" \"IPR000666\"";
		String iprAnno2 = "GOAT IPR000999";
		// Test
		DomainScoreCalculator.parseDomainAnnotationLine(iprAnno1);
		DomainScoreCalculator.parseDomainAnnotationLine(iprAnno2);
		assertTrue("Should have parsed two InterPro Domain Annotations.",
				!DomainScoreCalculator.getBlastResultAccessionsToInterproIds()
						.isEmpty());
		assertTrue("IPR000666 should have been annotated",
				DomainScoreCalculator.getBlastResultAccessionsToInterproIds()
						.get("SHEEP").contains("IPR000666"));
		assertTrue("IPR000999 should have been annotated",
				DomainScoreCalculator.getBlastResultAccessionsToInterproIds()
						.get("GOAT").contains("IPR000999"));
		// Setup
		getSettings().setComputeDomainSimilarityOn("pfam");
		String pfamAnno1 = "\"SHEEP\" \"PF12345\"";
		String pfamAnno2 = "GOAT PF54321";
		// Test
		DomainScoreCalculator.parseDomainAnnotationLine(pfamAnno1);
		DomainScoreCalculator.parseDomainAnnotationLine(pfamAnno2);
		assertTrue("Should have parsed two PFam Domain Annotations.",
				!DomainScoreCalculator.getBlastResultAccessionsToPfamIds()
						.isEmpty());
		assertTrue(
				"PF12345 should have been annotated",
				DomainScoreCalculator.getBlastResultAccessionsToPfamIds()
						.get("SHEEP").contains("PF12345"));
		assertTrue(
				"PF54321 should have been annotated",
				DomainScoreCalculator.getBlastResultAccessionsToPfamIds()
						.get("GOAT").contains("PF54321"));
	}
}
