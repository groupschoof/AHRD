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
import ahrd.model.Protein;
import nu.xom.ParsingException;

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

}
