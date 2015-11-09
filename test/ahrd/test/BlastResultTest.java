package ahrd.test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import static ahrd.controller.Settings.getSettings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import org.junit.Test;
import org.xml.sax.SAXException;

import ahrd.exception.MissingProteinException;
import ahrd.model.BlastResult;
import ahrd.model.Protein;

public class BlastResultTest {

	public BlastResultTest() {
		super();
	}

	@Test
	public void testPatternize() {
		BlastResult br = TestUtils.mockBlastResult("accession_5", 5.0, "description_5 Fly-Wing formation", 10, 20, 10,
				20, 200, 30.0, "trembl",
				new HashSet<String>(Arrays.asList("description", "5", "fly", "wing", "formation")));
		assertEquals("5descriptionflyformationwing", br.patternize());
	}

	@Test
	public void testTokenize() throws IOException {
		TestUtils.initTestSettings();
		BlastResult br = new BlastResult("accession_1", 1.0, "one tWo Three protein homolog putative", 10, 20, 10, 20,
				200, 30, "swissprot");
		br.tokenize();
		// test:
		assertEquals(3.0, br.getTokens().size(), 0.0);
		assertTrue(br.getTokens().contains("one"));
		assertTrue(br.getTokens().contains("two"));
		assertTrue(br.getTokens().contains("three"));
	}

	/**
	 * TODO: Remove overlapping tests with AhrdTest.testAhrdParsesBlast
	 * 
	 * @throws MissingProteinException
	 * @throws SAXException
	 * @throws IOException
	 */
	@Test
	public void testParseBlastResultsOfQueryProteins() throws MissingProteinException, SAXException, IOException {
		// Strangely the following line missing causes a NPE in Java 7:
		TestUtils.initTestSettings();
		Map<String, Protein> protDb = TestUtils.mockProteinDb();
		assertTrue(protDb.containsKey("gene:chr01.1056:mRNA:chr01.1056"));
		assertTrue(protDb.containsKey("gene:chr01.502:mRNA:chr01.502"));

		BlastResult.parseBlastResults(protDb, "swissprot");
		assertNotNull(protDb.get("gene:chr01.1056:mRNA:chr01.1056").getBlastResults().get("swissprot"));
		assertEquals(108, protDb.get("gene:chr01.1056:mRNA:chr01.1056").getBlastResults().get("swissprot").size());
		assertEquals(16, protDb.get("gene:chr01.502:mRNA:chr01.502").getBlastResults().get("swissprot").size());

		BlastResult br = (BlastResult) protDb.get("gene:chr01.1056:mRNA:chr01.1056").getBlastResults().get("swissprot")
				.get(0);
		assertEquals("sp|Q9SCZ4|FERON_ARATH", br.getAccession());
		assertEquals(0.0, br.getEValue(), 0.0);
		assertEquals(1095, br.getBitScore(), 0.0);
		assertEquals(30.0, br.getQueryStart(), 0.0);
		assertEquals(828.0, br.getQueryEnd(), 0.0);
		assertEquals("Receptor-like protein kinase FERONIA", br.getDescription());
		// While parsing the BlastResults for a Protein, the
		// frequencies of each Description-Line should be
		// measured. This is after the Description-Line
		// passes the blacklist-check and has it's "bad
		// tokens" filtered out.
		Protein p = protDb.get("gene:chr01.502:mRNA:chr01.502");
		// Test for maximum Bit-Score being saved:
		assertEquals(94.4, p.getDescriptionScoreCalculator().getMaxBitScore(), 0.0);
	}

	@Test
	public void testCompareBlastResultsBasedOnTheirEvalues() {
		BlastResult br1 = new BlastResult("accession_1", 1.0, "description_1", 10, 20, 10, 20, 200, 30, "swissprot");
		BlastResult br2 = new BlastResult("accession_2", 2.0, "description_2", 10, 20, 10, 20, 200, 30, "swissprot");
		BlastResult br3 = new BlastResult("accession_3", 2.0, "description_3", 10, 20, 10, 20, 200, 30, "swissprot");
		assertEquals(-1, br1.compareTo(br2));
		assertEquals(1, br2.compareTo(br1));
		assertEquals(0, br2.compareTo(br3));
	}

	@Test
	public void testFilterBestScoringBlastResults() {
		List<BlastResult> blastResults = new ArrayList<BlastResult>();
		blastResults.add(new BlastResult("accession_1", 1.0, "description_1", 10, 20, 10, 20, 200, 30, "swissprot"));
		blastResults.add(new BlastResult("accession_2", 2.0, "description_2", 10, 20, 10, 20, 200, 30, "swissprot"));
		blastResults.add(new BlastResult("accession_3", 3.0, "description_3", 10, 20, 10, 20, 200, 30, "swissprot"));
		BlastResult br = new BlastResult("accession_4", 4.0, "description_4", 10, 20, 10, 20, 200, 30, "swissprot");
		blastResults.add(br);
		blastResults.add(new BlastResult("accession_5", 5.0, "description_5", 10, 20, 10, 20, 200, 30, "swissprot"));
		blastResults = BlastResult.filterBestScoringBlastResults(blastResults, 4);
		assertEquals(4, blastResults.size());
		assertTrue(blastResults.contains(br));
	}

	@Test
	public void testGenerateHRDCandidateForProtein() throws IOException, MissingProteinException, SAXException {
		TestUtils.initTestSettings();
		getSettings().setWriteBestBlastHitsToOutput(true);
		Map<String, Protein> protDb = TestUtils.mockProteinDb();
		assertTrue(protDb.containsKey("gene:chr01.1056:mRNA:chr01.1056"));
		assertTrue(protDb.containsKey("gene:chr01.502:mRNA:chr01.502"));

		BlastResult.parseBlastResults(protDb, "swissprot");
		assertNotNull(protDb.get("gene:chr01.1056:mRNA:chr01.1056").getBlastResults().get("swissprot"));
		assertEquals(108, protDb.get("gene:chr01.1056:mRNA:chr01.1056").getBlastResults().get("swissprot").size());
		assertEquals(16, protDb.get("gene:chr01.502:mRNA:chr01.502").getBlastResults().get("swissprot").size());

		Protein p1 = protDb.get("gene:chr01.502:mRNA:chr01.502");
		assertNotNull(p1.getEvaluationScoreCalculator().getUnchangedBlastResults());
		assertTrue(!p1.getEvaluationScoreCalculator().getUnchangedBlastResults().isEmpty());
		assertEquals("sp|Q3EBC8|DCL2_ARATH",
				p1.getEvaluationScoreCalculator().getUnchangedBlastResults().get("swissprot").getAccession());

		Protein p2 = protDb.get("gene:chr01.1056:mRNA:chr01.1056");
		assertNotNull(p2.getEvaluationScoreCalculator().getUnchangedBlastResults());
		assertTrue(!p2.getEvaluationScoreCalculator().getUnchangedBlastResults().isEmpty());
		assertEquals("sp|Q9SCZ4|FERON_ARATH",
				p2.getEvaluationScoreCalculator().getUnchangedBlastResults().get("swissprot").getAccession());
	}
}
