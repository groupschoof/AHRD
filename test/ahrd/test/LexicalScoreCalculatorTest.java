package ahrd.test;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.junit.Before;
import org.junit.Test;

import ahrd.model.BlastResult;
import ahrd.model.Protein;

public class LexicalScoreCalculatorTest {

	@Before
	public void setUp() throws IOException {
		TestUtils.initTestSettings();
	}

	@Test
	public void testSumTokenScoresDividebByHighScore() {
		Protein p = TestUtils.mockProtein();
		List<BlastResult> brs = new ArrayList<BlastResult>();
		// Mocked BlastResult has Tokens one, two, three
		BlastResult br = TestUtils.mockBlastResult();
		brs.add(br);
		p.getBlastResults().put("swissprot", TestUtils.mockBlastResults());
		p.getTokenScoreCalculator().getTokenScores().put("one", 0.2);
		p.getTokenScoreCalculator().getTokenScores().put("two", 0.3);
		p.getTokenScoreCalculator().getTokenScores().put("three", 0.8);
		p.getTokenScoreCalculator().setTokenHighScore(0.8);
		// test: (0.2 + 0.3 + 0.8) / 0.8
		assertEquals(1.625, p.getLexicalScoreCalculator()
				.sumTokenScoresDividebByHighScore(br), 0.0);
	}

	@Test
	public void testCorrectionFactor() {
		Protein p = TestUtils.mockProtein();
		List<BlastResult> brs = new ArrayList<BlastResult>();
		// Mocked BlastResult has Tokens one, two, three
		BlastResult br = TestUtils.mockBlastResult();
		brs.add(br);
		p.getBlastResults().put("swissprot", TestUtils.mockBlastResults());
		p.getTokenScoreCalculator().getTokenScores().put("one", 0.222);
		p.getTokenScoreCalculator().getTokenScores().put("two", 0.333);
		p.getTokenScoreCalculator().getTokenScores().put("three", 0.888);
		p.getTokenScoreCalculator().setTokenHighScore(0.888);

		assertEquals((3.0 / 1.0), p.getLexicalScoreCalculator()
				.correctionFactor(br), 0.0);
	}

	@Test
	public void testLexicalScore() {
		Protein p = TestUtils.mockProtein();
		List<BlastResult> brs = new ArrayList<BlastResult>();
		BlastResult br = TestUtils.mockBlastResult();
		brs.add(br);
		// GeneOntologyScore is always 1.11
		p.setLexicalScoreCalculator(new TestUtils.LexicalScoreCalculatorFixedGoScoreMock(
				p));
		// Mocked BlastResult has Tokens one, two, three
		p.getBlastResults().put("swissprot", TestUtils.mockBlastResults());
		// Note: Token-Scores are NOT filtered in this test!
		p.getTokenScoreCalculator().getTokenScores().put("one", 0.2);
		p.getTokenScoreCalculator().getTokenScores().put("two", 0.3);
		p.getTokenScoreCalculator().getTokenScores().put("three", 0.8);
		p.getTokenScoreCalculator().setTokenHighScore(0.8);

		// (((0.2 + 0.3 + 0.8) / 0.8) / (3.0/1.0))
		assertEquals(0.5416667, p.getLexicalScoreCalculator().lexicalScore(br),
				0.0000001);
	}
}
