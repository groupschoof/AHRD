package ahrd.test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.io.IOException;

import org.junit.Before;
import org.junit.Test;
import org.xml.sax.SAXException;

import ahrd.controller.AHRD;
import ahrd.exception.MissingAccessionException;
import ahrd.exception.MissingProteinException;
import ahrd.model.BlastResult;
import ahrd.model.InterproResult;
import ahrd.model.Protein;

public class AhrdTest {

	private AHRD ahrd;

	public AhrdTest() {
		super();
	}

	@Before
	public void setUp() throws IOException {
		ahrd = new AHRD("./test/resources/ahrd_input.yml");
	}

	@Test
	public void testAhrdInitializesProteins() throws IOException,
			MissingAccessionException {
		ahrd.initializeProteins();
		assertNotNull(ahrd.getProteins());
		assertTrue(ahrd.getProteins().containsKey(
				"gene:chr01.502:mRNA:chr01.502"));
		assertTrue(ahrd.getProteins().containsKey(
				"gene:chr01.1056:mRNA:chr01.1056"));
	}

	@Test
	public void testAhrdParsesBlast() throws IOException,
			MissingProteinException, SAXException, MissingAccessionException {
		// We need the test-protein-Database in memory:
		ahrd.setProteins(TestUtils.mockProteinDb());
		ahrd.parseBlastResults();
		for (String iterProtAcc : ahrd.getProteins().keySet()) {
			Protein protein = ahrd.getProteins().get(iterProtAcc);
			assertTrue(protein.getBlastResults().size() > 0);
			assertTrue(protein.getBlastResults().containsKey("swissprot"));
			assertTrue(protein.getBlastResults().containsKey("tair"));
			assertTrue(protein.getBlastResults().containsKey("trembl"));
		}
		// test number of tokens:
		Protein p = ahrd.getProteins().get("gene:chr01.502:mRNA:chr01.502");
		BlastResult br = p.getBlastResults().get("swissprot").get(0);
		assertEquals(2.0, br.getTokens().size(), 0.0);
		// test measurement of cumulative token-scores was triggered:
		assertTrue(p.getTokenScoreCalculator().getCumulativeTokenBitScores()
				.containsKey("dicer"));
		assertTrue(p.getTokenScoreCalculator()
				.getCumulativeTokenBlastDatabaseScores().containsKey("dicer"));
		assertTrue(p.getTokenScoreCalculator()
				.getCumulativeTokenOverlapScores().containsKey("dicer"));
		// test measurement of total token-scores was triggered:
		assertTrue(p.getTokenScoreCalculator().getTotalTokenBitScore() > 0.0);
		assertTrue(p.getTokenScoreCalculator()
				.getTotalTokenBlastDatabaseScore() > 0.0);
		assertTrue(p.getTokenScoreCalculator().getTotalTokenOverlapScore() > 0.0);
	}

	@Test
	public void testParseInterproResults() throws Exception {
		ahrd.setProteins(TestUtils.mockProteinDb());
		// TODO: The Interpro-Database should be mocked:
		InterproResult.initialiseInterproDb();
		ahrd.parseInterproResult();
		for (String iterProtAcc : ahrd.getProteins().keySet()) {
			Protein protein = ahrd.getProteins().get(iterProtAcc);
			assertTrue(protein.getInterproResults().size() > 0);
		}
	}

	@Test
	public void testParseGeneOntologyResults() throws Exception {
		ahrd.setProteins(TestUtils.mockProteinDb());
		ahrd.parseGeneOntologyResult();
		assertEquals(1, ahrd.getProteins().get("gene:chr01.502:mRNA:chr01.502")
				.getGoResults().size());
		assertEquals(2, ahrd.getProteins().get(
				"gene:chr01.1056:mRNA:chr01.1056").getGoResults().size());
	}

	// @Test
	// public void testScoresInCompleteRun() throws IOException,
	// MissingAccessionException, MissingProteinException, SAXException,
	// ParsingException {
	// System.out.println("Method: testScoresInCompleteRun()");
	// // Setup a AHRD-Run for a single Protein-Sequence
	// AHRD a = new AHRD("./test/resources/scores_test/ahrd_input.yml");
	// a.setup(false);
	// // Test measurements:
	// Protein p = a.getProteins().get("Solyc11g010630.2.1");
	// assertEquals(1.0, p.getDescriptionScoreCalculator()
	// .getDescriptionLineFrequencies().get("AT5g65480/K19O4_1"), 0.0);
	// assertEquals(1.0,
	// p.getDescriptionScoreCalculator().getMaxDescriptionLineFrequency(), 0.0);
	// // Find highest scoring Description-Line and test for correct scores:
	// a.filterBestScoringBlastResults(p);
	// assertEquals(1.0, p.getBlastResults().get("trembl").size(), 0.0);
	// // a.tokenizeBlastResultDescriptionLines(p);
	//		
	// // Only a single BlastHit in trembl:
	// assertEquals(2.0,
	// p.getBlastResults().get("trembl").get(0).getTokens().size(), 0.0);
	// // Only trembl has processable BlastHits:
	// p.getTokenScoreCalculator().assignTokenScores();
	// assertEquals(2.0, p.getTokenScoreCalculator().getTokenScores().size(),
	// 0.0);
	// assertEquals(1.0,
	// p.getTokenScoreCalculator().getTokenScores().get("at5g65480"), 0.0);
	// assertEquals(1.0,
	// p.getTokenScoreCalculator().getTokenScores().get("k1904_1"), 0.0);
	// // prot.getTokenScoreCalculator().filterTokenScores();
	// }

}
