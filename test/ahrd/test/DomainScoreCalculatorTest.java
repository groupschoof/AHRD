package ahrd.test;

import static junit.framework.Assert.*;

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

import org.junit.Before;
import org.junit.Test;

import ahrd.model.BlastResult;
import ahrd.model.DomainScoreCalculator;
import ahrd.model.InterproResult;
import ahrd.model.Protein;

public class DomainScoreCalculatorTest {

	@Before
	public void initialiseInterproDb() throws IOException {
		TestUtils.initTestSettings();
	}

	@Test
	public void testInitializeBlastResultAccessionsToInterproIds()
			throws IOException {
		DomainScoreCalculator.initializeBlastResultAccessionsToInterproIds();
		assertTrue("IPR012610 should be assigned to DCL2_ARATH",
				DomainScoreCalculator.getBlastResultAccessionsToInterproIds()
						.get("DCL2_ARATH").contains("IPR012610"));
		assertTrue("IPR012610 should be assigned to DCL2A_ORYSJ",
				DomainScoreCalculator.getBlastResultAccessionsToInterproIds()
						.get("DCL2A_ORYSJ").contains("IPR012610"));
		assertTrue("IPR020139 should be assigned to DCL2_ARATH",
				DomainScoreCalculator.getBlastResultAccessionsToInterproIds()
						.get("DCL2_ARATH").contains("IPR020139"));
	}

	@Test
	public void testConstructVectorSpaceModel() {
		// Protein to be tested:
		Protein prot = TestUtils.mockProtein();
		// Interpro-Domains, it has been annotated with:
		prot.getInterproResults().add(
				new InterproResult("IPR000001", "first", "Domain"));
		prot.getInterproResults().add(
				new InterproResult("IPR000002", "second", "Domain"));
		prot.getInterproResults().add(
				new InterproResult("IPR000003", "third", "Domain"));
		// BlastResult found in Uniprot/Swissprot for above Protein,
		// with accessions named 'accession_X' X = 1,2,...,5
		prot.getBlastResults().put("swissprot", TestUtils.mockBlastResults());
		// A single hit in TAIR:
		List<BlastResult> tairHits = new ArrayList<BlastResult>();
		tairHits.add(TestUtils.mockBlastResult("accession_6", 0.001,
				"description six", 1, 20, 100.0, "tair", new HashSet<String>()));
		prot.getBlastResults().put("tair", tairHits);
		// Setup memory database 'BlastResultAccessionsToInterproIds'
		String interpro1 = "IPR000001";
		String interpro2 = "IPR000002";
		String interpro3 = "IPR000003";
		String interpro4 = "IPR000004";
		String interpro5 = "IPR000005";
		String interpro6 = "IPR000006";
		Map<String, Set<String>> brAccs2iprIds = new HashMap<String, Set<String>>();
		brAccs2iprIds.put(
				"accession_1",
				new HashSet<String>(Arrays.asList(new String[] { interpro1,
						interpro2 })));
		brAccs2iprIds.put(
				"accession_2",
				new HashSet<String>(Arrays.asList(new String[] { interpro1,
						interpro3 })));
		brAccs2iprIds.put("accession_3",
				new HashSet<String>(Arrays.asList(new String[] { interpro4 })));
		brAccs2iprIds.put("accession_4",
				new HashSet<String>(Arrays.asList(new String[] { interpro5 })));
		// accession_5 has no InterproIDs assigned!
		brAccs2iprIds.put(
				"accession_6",
				new HashSet<String>(Arrays.asList(new String[] { interpro1,
						interpro2, interpro6 })));
		DomainScoreCalculator
				.setBlastResultAccessionsToInterproIds(brAccs2iprIds);

		// Test the construction of the vector space model:
		SortedSet<String> vectorSpaceModelToTest = DomainScoreCalculator
				.constructVectorSpaceModel(prot);
		SortedSet<String> vectorSpaceModelExpected = new TreeSet<String>(
				Arrays.asList(new String[] { interpro1, interpro2, interpro3,
						interpro4, interpro5, interpro6 }));
		assertEquals(
				"Construction of vector space model should return the following dimension in their alphabetical order: IPR00000X, X=1,2,...,6",
				vectorSpaceModelExpected, vectorSpaceModelToTest);

	}
}
