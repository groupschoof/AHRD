package ahrd.test;

import static org.junit.Assert.fail;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static ahrd.controller.Settings.getSettings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.junit.Before;
import org.junit.Test;

import ahrd.model.BlastResult;
import ahrd.model.Protein;
import ahrd.model.TokenScoreCalculator;

public class TokenScoreCalculatorTest {

	private BlastResult br1, br2, br3;

	public TokenScoreCalculatorTest() {
		super();
	}

	@Before
	public void setup() throws IOException {
		TestUtils.initTestSettings();
		br1 = new BlastResult("accession_1", 1.0, "description_1", 40, 99, 40,
				99, 200, 69.96, "swissprot");
		br1.getTokens().add("token_one");
		// ensure double tokens have no effect:
		br1.getTokens().add("token_one");
		br2 = new BlastResult("accession_2", 2.0, "description_2", 45, 69, 45,
				69, 200, 45.54, "tair");
		br2.getTokens().add("token_one");
		br2.getTokens().add("token_two");
		br3 = new BlastResult("accession_3", 2.0, "description_3", 35, 125, 35,
				125, 200, 88.0, "trembl");
		br3.getTokens().add("token_two");
		br3.getTokens().add("token_three");
	}

	@Test
	public void testOverlapScore() {
		double subjectEnd = 100;
		double subjectStart = 10;
		double subjectLength = 200;
		double queryEnd = 190;
		double queryStart = 110;
		double queryLength = 200;
		// (100 - 10 + 190 - 110 + 2) / 400 = 0.43
		assertEquals(0.43,
				TokenScoreCalculator.overlapScore(queryStart, queryEnd,
						queryLength, subjectStart, subjectEnd, subjectLength),
				0.0000001);
	}

	@Test
	public void testMeasureTotalScores() throws IOException {
		// Init the test-settings:
		TestUtils.initTestSettings();
		Protein p = TestUtils.mockProtein(); // sequence-length = 200
		p.getTokenScoreCalculator().measureTotalScores(br1);
		p.getTokenScoreCalculator().measureTotalScores(br2);
		p.getTokenScoreCalculator().measureTotalScores(br3);
		TokenScoreCalculator tsc = p.getTokenScoreCalculator();

		// test
		assertEquals(203.5, tsc.getTotalTokenBitScore(), 0.0);
		assertEquals(160.0, tsc.getTotalTokenBlastDatabaseScore(), 0.0);
		assertEquals(0.88, tsc.getTotalTokenOverlapScore(), 0.0);
	}

	@Test
	public void testMeasureCumulativeScores() throws IOException {
		// Init the test-settings:
		TestUtils.initTestSettings();
		Protein p = TestUtils.mockProtein(); // sequence-length = 200
		p.getTokenScoreCalculator().measureCumulativeScores(br1);
		p.getTokenScoreCalculator().measureCumulativeScores(br2);
		p.getTokenScoreCalculator().measureCumulativeScores(br3);
		TokenScoreCalculator tsc = p.getTokenScoreCalculator();

		// We have three tokens:
		assertEquals(3, tsc.getCumulativeTokenBitScores().size());
		assertEquals(3, tsc.getCumulativeTokenOverlapScores().size());
		assertEquals(3, tsc.getCumulativeTokenBlastDatabaseScores().size());
		// test cum.BitScores
		assertEquals(115.5, tsc.getCumulativeTokenBitScores().get("token_one"),
				0);
		assertEquals(133.54,
				tsc.getCumulativeTokenBitScores().get("token_two"), 0);
		assertEquals(88.0,
				tsc.getCumulativeTokenBitScores().get("token_three"), 0);
		// test cum.BlastDatabaseScore
		assertEquals(150,
				tsc.getCumulativeTokenBlastDatabaseScores().get("token_one"), 0);
		assertEquals(60,
				tsc.getCumulativeTokenBlastDatabaseScores().get("token_two"), 0);
		assertEquals(10,
				tsc.getCumulativeTokenBlastDatabaseScores().get("token_three"),
				0);
		// test cum.OverlapScores:
		assertEquals(0.425,
				tsc.getCumulativeTokenOverlapScores().get("token_one"), 0);
		assertEquals(0.5800000000000001, tsc.getCumulativeTokenOverlapScores()
				.get("token_two"), 0);
		assertEquals(0.455,
				tsc.getCumulativeTokenOverlapScores().get("token_three"), 0);
	}

	@Test
	public void testSumOfAllTokenScores() {
		Protein p = TestUtils.mockProtein();
		BlastResult one = new BlastResult("accession_1", 1.0, "first", 40, 99,
				40, 99, 200, 69.96, "swissprot");
		BlastResult two = new BlastResult("accession_2", 2.0, "first second",
				45, 69, 45, 69, 200, 45.54, "tair");
		one.getTokens().add("first");
		two.getTokens().add("first");
		two.getTokens().add("second");
		p.getTokenScoreCalculator().getTokenScores().put("first", 0.25);
		p.getTokenScoreCalculator().getTokenScores().put("second", 0.50);
		assertEquals(0.25,
				p.getTokenScoreCalculator().sumOfAllTokenScores(one), 0);
		assertEquals(0.75,
				p.getTokenScoreCalculator().sumOfAllTokenScores(two), 0);
	}

	@Test
	public void testAssignTokenScores() {
		// Mock Protein
		Protein p = TestUtils.mockProtein();
		// Mock BlastResults for Protein p,
		// where Tokens are: one, two, three
		List<BlastResult> mockedBlastResults = TestUtils
				.mockBlastResultsWithTokens();
		// Set BlasteRsults for two BlastDatabases
		Map<String, List<BlastResult>> blastResults = new HashMap<String, List<BlastResult>>();
		blastResults.put("swissprot", mockedBlastResults);
		blastResults.put("tair", mockedBlastResults);
		p.setBlastResults(blastResults);
		// Mock "Running-Sums" and "Totals", which are set by method
		// "TokenScoreCalculator.addToken(..)"
		TestUtils.mockCumulativeTokenScores(p, "one", 2);
		TestUtils.mockCumulativeTokenScores(p, "two", 5);
		TestUtils.mockCumulativeTokenScores(p, "three", 10);
		p.getTokenScoreCalculator().setTotalTokenBitScore(250);
		p.getTokenScoreCalculator().setTotalTokenBlastDatabaseScore(300);
		p.getTokenScoreCalculator().setTotalTokenOverlapScore(3.75);

		// Call method to test:
		p.getTokenScoreCalculator().assignTokenScores();
		// Assert expectations
		assertTrue(p.getTokenScoreCalculator().getTokenScores()
				.containsKey("one"));
		assertTrue(p.getTokenScoreCalculator().getTokenScores()
				.containsKey("two"));
		assertTrue(p.getTokenScoreCalculator().getTokenScores()
				.containsKey("three"));
		assertTrue(p.getTokenScoreCalculator().getTokenHighScore() > 0);
		assertEquals(0.22666666666666668, p.getTokenScoreCalculator()
				.getTokenHighScore(), 0.0000001);
	}

	@Test
	public void testFilterTokenScores() {
		Protein p = TestUtils.mockProtein();
		TokenScoreCalculator tsc = p.getTokenScoreCalculator();
		tsc.setTokenHighScore(0.666);
		tsc.getTokenScores().put("sheep", 0.222);
		tsc.getTokenScores().put("goat", 0.444);
		tsc.getTokenScores().put("ram", 0.111);
		tsc.getTokenScores().put("batsheep", 0.555);
		tsc.filterTokenScores();
		// 0.666 / 2 = 0.333
		assertEquals(-0.11100000000000002, tsc.getTokenScores().get("sheep"), 0);
		assertEquals((0.444), tsc.getTokenScores().get("goat"), 0);
		assertEquals(-0.22200000000000003, tsc.getTokenScores().get("ram"), 0);
		assertEquals((0.555), tsc.getTokenScores().get("batsheep"), 0);
	}

	@Test
	public void testDescriptionLineSummedTokenScore() {
		Protein p = TestUtils.mockProtein();
		List<BlastResult> brs = new ArrayList<BlastResult>();
		// Mocked BlastResult has Tokens one, two, three
		brs.add(TestUtils.mockBlastResult());
		p.getBlastResults().put("swissprot", brs);
		// Add some arbitrary scores:
		p.getTokenScoreCalculator().getTokenScores().put("one", 0.222);
		p.getTokenScoreCalculator().getTokenScores().put("two", 0.333);
		p.getTokenScoreCalculator().getTokenScores().put("three", 0.445);

		assertEquals(1.0, p.getTokenScoreCalculator()
				.descriptionLineSummedTokenScore(brs.get(0)), 0.0);
	}

	@Test
	public void testValidationOfTokenScoreWeights() {
		Protein p = TestUtils.mockProtein();
		TokenScoreCalculator tsc = p.getTokenScoreCalculator();
		String token = "foo";
		// Init cumulative scores:
		tsc.getCumulativeTokenBitScores().put(token, 0.5);
		tsc.getCumulativeTokenBlastDatabaseScores().put(token, 0.5);
		tsc.getCumulativeTokenOverlapScores().put(token, 0.5);

		getSettings().setDescriptionTokenScoreBitScoreWeight(0.5);
		getSettings().setDescriptionTokenScoreDatabaseScoreWeight(0.5);
		getSettings().setDescriptionTokenScoreOverlapScoreWeight(0.0011);
		try {
			tsc.tokenScore("foo", "swissprot");
			fail("Validation of the three weights in the formula Token-Score failed. Their sum should be >= 0.999 and <= 1.001");
		} catch (IllegalArgumentException expectedException) {
		}

		getSettings().setDescriptionTokenScoreBitScoreWeight(0.5);
		getSettings().setDescriptionTokenScoreDatabaseScoreWeight(0.3);
		getSettings().setDescriptionTokenScoreOverlapScoreWeight(0.198);
		try {
			tsc.tokenScore(token, "swissprot");
			fail("Validation of the three weights in the formula Token-Score failed. Their sum should be >= 0.999 and <= 1.001");
		} catch (IllegalArgumentException expectedException) {
		}

		// Mock:
		tsc.getCumulativeTokenBitScores().put("foo", 0.5);
		tsc.getCumulativeTokenBlastDatabaseScores().put("foo", 0.6);
		tsc.getCumulativeTokenOverlapScores().put("foo", 0.7);

		getSettings().setDescriptionTokenScoreBitScoreWeight(0.5);
		getSettings().setDescriptionTokenScoreDatabaseScoreWeight(0.5);
		getSettings().setDescriptionTokenScoreOverlapScoreWeight(0.001);
		try {
			tsc.tokenScore(token, "swissprot");
		} catch (IllegalArgumentException expectedException) {
			fail("Validation of the three weights in the formula Token-Score failed. It is too restrictive, a delta of 0.001 has to be excepted.");
		}

		getSettings().setDescriptionTokenScoreBitScoreWeight(0.5);
		getSettings().setDescriptionTokenScoreDatabaseScoreWeight(0.3);
		getSettings().setDescriptionTokenScoreOverlapScoreWeight(0.199);
		try {
			tsc.tokenScore(token, "swissprot");
		} catch (IllegalArgumentException expectedException) {
			fail("Validation of the three weights in the formula Token-Score failed. It is too restrictive, a delta of 0.001 has to be excepted.");
		}
	}
}
