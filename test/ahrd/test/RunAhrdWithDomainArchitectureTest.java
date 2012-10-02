package ahrd.test;

import static ahrd.controller.Settings.getSettings;
import static ahrd.controller.Utils.roundEachToNDecimalPlaces;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.Set;

import org.junit.Before;
import org.junit.Test;

import ahrd.controller.AHRD;
import ahrd.model.BlastResult;
import ahrd.model.DomainScoreCalculator;
import ahrd.model.InterproResult;
import ahrd.model.Protein;

public class RunAhrdWithDomainArchitectureTest {

	private AHRD ahrd;

	@Before
	public void setUp() throws IOException {
		ahrd = new AHRD("./test/resources/ahrd_dom_arch_sim_input.yml");
	}

	/**
	 * Check each step of AHRD executed with options set to take domain
	 * architecture similarity (DAS) into account. DAS is computed on InterPro
	 * domains.
	 * 
	 * @throws Exception
	 */
	@Test
	public void testRunInterpro() throws Exception {

		// In this TEST BLOCK ignore the domain architecture similarity score as
		// a summand for the token score:
		getSettings().setTokenScoreDomainSimilarityWeight(0.0);

		// Tell User about this ToDo:
		System.err
				.println("WARNING: There is still no test implemented to validate AHRD's behaviour when used with a non zero weight for 'TokenScoreDomainSimilarityWeight'."
						+ "\nPlease implement such a test and delete this warning message.");

		// Check correct initialization:
		ahrd.setup(false);

		// Verify the expected InterPro domains have been loaded:
		assertNotNull(
				"InterPro Entry 'IPR000001' was not correctly loaded from the interpro database!",
				InterproResult.getInterproDb().get("IPR000001"));
		assertNotNull(
				"InterPro Entry 'IPR000003' was not correctly loaded from the interpro database!",
				InterproResult.getInterproDb().get("IPR000003"));
		assertNotNull(
				"InterPro Entry 'IPR000006' was not correctly loaded from the interpro database!",
				InterproResult.getInterproDb().get("IPR000006"));
		assertNotNull(
				"InterPro Entry 'IPR000535' was not correctly loaded from the interpro database!",
				InterproResult.getInterproDb().get("IPR000535"));
		assertNotNull(
				"InterPro Entry 'IPR000536' was not correctly loaded from the interpro database!",
				InterproResult.getInterproDb().get("IPR000536"));
		assertNotNull(
				"InterPro Entry 'IPR001723' was not correctly loaded from the interpro database!",
				InterproResult.getInterproDb().get("IPR001723"));
		assertNotNull(
				"InterPro Entry 'IPR001732' was not correctly loaded from the interpro database!",
				InterproResult.getInterproDb().get("IPR001732"));
		assertNotNull(
				"InterPro Entry 'IPR001733' was not correctly loaded from the interpro database!",
				InterproResult.getInterproDb().get("IPR001733"));
		assertNotNull(
				"InterPro Entry 'IPR003019' was not correctly loaded from the interpro database!",
				InterproResult.getInterproDb().get("IPR003019"));
		assertNotNull(
				"InterPro Entry 'IPR008946' was not correctly loaded from the interpro database!",
				InterproResult.getInterproDb().get("IPR008946"));
		assertNotNull(
				"InterPro Entry 'IPR013806' was not correctly loaded from the interpro database!",
				InterproResult.getInterproDb().get("IPR013806"));
		assertNotNull(
				"InterPro Entry 'IPR016040' was not correctly loaded from the interpro database!",
				InterproResult.getInterproDb().get("IPR016040"));
		assertNotNull(
				"InterPro Entry 'IPR019756' was not correctly loaded from the interpro database!",
				InterproResult.getInterproDb().get("IPR019756"));

		// Domain Weights Database has to have been parsed:
		InterproResult ipr = InterproResult.getInterproDb().get("IPR000001");
		Double dw = ipr.getDomainWeight();
		assertNotNull(
				"InterPro Entry 'IPR000001' has no Domain Weight assigned!", dw);
		assertEquals(0.189433136086726, dw, 0.0);

		// Test AHRD:
		ahrd.assignHumanReadableDescriptions();

		// With InterPro domain annotations:
		Protein p1 = ahrd.getProteins().get("gene:chr01.502:mRNA:chr01.502");
		BlastResult bestBr1 = p1.getDescriptionScoreCalculator()
				.getHighestScoringBlastResult();
		Double descScore1 = bestBr1.getDescriptionScore();
		assertNotNull(
				"Description Score of Blast Hit 'sp|Q3EBC8|DCL2_ARATH' should not be NULL!",
				descScore1);
		assertNotNull(
				"Protein 'gene:chr01.502:mRNA:chr01.502' should have a vector in domain architecture space.",
				p1.getDomainWeights());
		assertNotNull(
				"BlastResult 'sp|Q3EBC8|DCL2_ARATH' should have a vector in domain architecture space.",
				bestBr1.getDomainWeights());
		assertNotNull(
				"BlastResult 'sp|Q3EBC8|DCL2_ARATH' should have a computed Domain Architecture Similarity Score.",
				bestBr1.getDomainSimilarityScore());
		// But no domains are shared with the best BlastHit, so:
		assertEquals(0.0, bestBr1.getDomainSimilarityScore(), 0.0);
		// Expect Description Score to be equal to score calculated as without
		// domain annotations:
		assertEquals(2.947, bestBr1.getDescriptionScore(), 0.001);

		Protein p2 = ahrd.getProteins().get("Solyc11g030630.1.1");
		assertNotNull(
				"Protein 'Solyc11g030630.1.1' should have Domain Annotations.",
				p2.getInterproResults());
		assertEquals(6, p2.getInterproResults().size());
		BlastResult bestBr2 = p2.getDescriptionScoreCalculator()
				.getHighestScoringBlastResult();
		Set<String> bestBr2IprAnnos = DomainScoreCalculator
				.getBlastResultAccessionsToInterproIds().get(
						bestBr2.getAccession());
		assertNotNull(
				"BlastResult '"
						+ bestBr2.getAccession()
						+ "' should have been assigned InterPro Annotations from Uniprot.",
				bestBr2IprAnnos);
		assertTrue("BlastResult '" + bestBr2.getAccession()
				+ "' should have more than six Domain Annotations.",
				bestBr2IprAnnos.size() > 6);
		assertTrue(
				"BlastResult '"
						+ bestBr2.getAccession()
						+ "' and Query Protein '' should share at least a single Domain Annotation.",
				bestBr2IprAnnos.contains(((InterproResult) p2
						.getInterproResults().toArray()[0]).getId()));
		Double descScore2 = bestBr2.getDescriptionScore();
		assertNotNull(
				"Description Score of Blast Hit '" + bestBr2.getAccession()
						+ "' should not be NULL!", descScore2);
		assertNotNull(
				"Protein 'Solyc11g030630.1.1' should have a vector in domain architecture space.",
				p2.getDomainWeights());
		// 'p2' has the following Domain Annotations:
		// IPR000999, IPR001159, IPR001650, IPR003100, IPR005034, IPR011545
		List<Double> expectedDomainWeights = roundEachToNDecimalPlaces(
				Arrays.asList(new Double[] { 5.05047155114809,
						0.230164198497054, 0.0244178837367414,
						0.33175852306605, 1.6165065434774, 0.173963313549586 }),
				4);
		List<Double> actualDomainWeights = roundEachToNDecimalPlaces(
				p2.getDomainWeights(), 4);
		assertTrue(
				"The vector in domain architecture space of Protein '"
						+ p2.getAccession()
						+ "' should contain the domain weights for its annotated domains.",
				actualDomainWeights.containsAll(expectedDomainWeights));
		// bestBr2 has the following Domain Annotations:
		// IPR000999, IPR001159, IPR001650, IPR003100, IPR005034, IPR011545,
		// IPR014001
		assertNotNull("BlastResult '" + bestBr2.getAccession()
				+ "' should have a vector in domain architecture space.",
				bestBr2.getDomainWeights());
		assertTrue(
				"The vector in domain architecture space of BlastResult '"
						+ bestBr2.getAccession()
						+ "' should contain the domain weights of those domains it has been annotated with",
				bestBr2.getDomainWeights().containsAll(
						Arrays.asList(new Double[] { 5.05047155114809,
								0.230164198497054, 0.0244178837367414,
								0.33175852306605, 1.6165065434774,
								0.173963313549586, 0.0 })));
		// Hence the Vector Space Model should be the sorted List of bestBrs'
		// Domain Annotations!
		assertNotNull(
				"The Protein's DomainScoreCalculator should have been assigned a Vector Space Model of all distinct annotated domains.",
				p2.getDomainScoreCalculator().getVectorSpaceModel());
		assertTrue(
				"The generated vector space model should contain all domains the protein has been annotated with.",
				p2.getDomainScoreCalculator()
						.getVectorSpaceModel()
						.containsAll(
								Arrays.asList(new String[] { "IPR000999",
										"IPR001159", "IPR001650", "IPR003100",
										"IPR005034", "IPR011545", "IPR014001" })));
		// As the only domain bestBr2 has annotated and p2 hasn't and this
		// domain does not appear in the domain weights db, it should receive a
		// weight of 0.0:
		assertEquals(0.0, InterproResult.getInterproDb().get("IPR014001")
				.getDomainWeight(), 0.0);
		// ... and hence the Domain Architecture Similarity Score should be
		// maximal:
		assertNotNull(
				"BlastResult '"
						+ bestBr2.getAccession()
						+ "' should have a computed Domain Architecture Similarity Score.",
				bestBr2.getDomainSimilarityScore());
		assertEquals(0.0, bestBr2.getDomainSimilarityScore(), 1.0);
		// The correctness of score computation is checked elsewhere, here just
		// assure that sharing some conserved domains increases the AHRD-Score
		// above the "original" one without considering shared domain
		// architecture:
		assertTrue(
				"Sharing some domain architecture should increase AHRD's final description score above the score not taking domain architecture into account.",
				2.947 < bestBr2.getDescriptionScore());
	}
}