package ahrd.test;

import static ahrd.controller.Settings.getSettings;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;

import org.junit.Before;
import org.junit.Test;

import ahrd.model.BlastResult;
import ahrd.model.DescriptionScoreCalculator;
import ahrd.model.Protein;

public class DescriptionScoreCalculatorTest {

	public DescriptionScoreCalculatorTest() {
		super();
	}

	@Before
	public void setUp() throws IOException {
		TestUtils.initTestSettings();
	}

	@Test
	public void testPatternFactor() {
		Protein p = TestUtils.mockProtein();
		// br has description 'one two three'
		BlastResult br = TestUtils.mockBlastResult();
		DescriptionScoreCalculator dsc = p.getDescriptionScoreCalculator();
		dsc.getDescLinePatternFrequencies().put("sheepase", 5);
		dsc.getDescLinePatternFrequencies().put("onethreetwo", 10);
		dsc.setMaxDescriptionLineFrequency(10);
		// test:
		assertEquals(0.6, dsc.patternFactor(br), 0.0);
	}

	@Test
	public void testRelativeBlastScore() {
		// Mock test-data
		Protein p = TestUtils.mockProtein();
		// br has bit-score 30.0
		BlastResult br = TestUtils.mockBlastResult();
		p.getDescriptionScoreCalculator().setMaxBitScore(60.0);
		// test: (0.2 * 30.0 / 60.0)
		assertEquals(0.1,
				p.getDescriptionScoreCalculator().relativeBlastScore(br), 0.0);
	}

	@Test
	public void testDescriptionScore() {
		Protein p = new Protein("sweet_sheep_protein", 200);
		// and mock LexicalScoreCalculator
		p.setLexicalScoreCalculator(new TestUtils.LexicalScoreCalculatorMock(p));
		BlastResult br = TestUtils.mockBlastResult("accession", 1.0,
				"goat sheep wool", 10, 20, 10, 20, 200, 30.0, "swissprot",
				new HashSet<String>(Arrays.asList("goat", "sheep", "wool")));
		List<BlastResult> brs = new ArrayList<BlastResult>();
		brs.add(br);
		p.getBlastResults().put("swissprot", brs);
		p.getDescriptionScoreCalculator().setMaxBitScore(30);
		p.getDescriptionScoreCalculator().getDescLinePatternFrequencies()
				.put("goatsheepwool", 10);
		p.getDescriptionScoreCalculator().setMaxDescriptionLineFrequency(10);

		// Token-Scores are not needed, as the lexical score is mocked!
		// DescriptionScore(1.5) := mockedLexicalScore(0.70) + 0.6 *
		// PatternFactor(10/10) + 0.2 * BitScore(30/30)
		p.getDescriptionScoreCalculator().calcDescriptionScore(br);
		assertEquals(1.5, br.getDescriptionScore(), 0.0);

		// Test Description-Score with Interpro annotations, that is including a
		// non zero domain similarity score:
		getSettings().setDescriptionScoreDomainSimilarityWeight(0.5);
		br.setDomainSimilarityScore(0.9);
		p.getDescriptionScoreCalculator().setMaxDomainSimilarityScore(1.0);
		p.getDescriptionScoreCalculator().calcDescriptionScore(br);
		// DescriptionScore(1.95) := above_DescriptionScore(1.5) + 0.5 * 0.9/1.0
		assertEquals(
				"Having a domain similarity score set for the BlastResult its description score should include it.",
				1.95, br.getDescriptionScore(), 0.0);
	}

	@Test
	public void testFindHighestScoringBlastResult() {
		Protein p = TestUtils.mockProtein();
		// Sprot
		p.getBlastResults().put("swissprot",
				TestUtils.mockBlastResultsForDescCalcTest());
		// trEMBL
		p.getBlastResults().put(
				"trembl",
				Arrays.asList(TestUtils.mockBlastResult(
						"accession_5",
						5.0,
						"description_5 Fly-Wing formation",
						10,
						20,
						10,
						20,
						200,
						30.0,
						"trembl",
						new HashSet<String>(Arrays.asList("description", "5",
								"fly", "wing", "formation")))));
		p.setLexicalScoreCalculator(new TestUtils.LexicalScoreCalculatorMock(p));
		p.getDescriptionScoreCalculator().setMaxBitScore(30.0);
		p.getDescriptionScoreCalculator().setDescLinePatternFrequencies(
				TestUtils.mockDescriptionLineFrequenciesForDescCalcTest());
		p.getDescriptionScoreCalculator().setMaxDescriptionLineFrequency(5);
		// Token-Scores are not needed, as the lexical score is mocked!
		// test
		p.getDescriptionScoreCalculator().findHighestScoringBlastResult();
		// 0.7 (mocked) + 0.4 * 30/30 + 0.6 * 5/5
		assertEquals(1.7000000000000002, p.getDescriptionScoreCalculator()
				.getDescriptionHighScore(), 0.0);
		assertEquals("description_5 Fly-Wing formation", p
				.getDescriptionScoreCalculator().getHighestScoringBlastResult()
				.getDescription());
	}

	@Test
	public void testDomainSimilarityScore() {
		Protein p = TestUtils.mockProtein();
		BlastResult br = new BlastResult("accession_1", 0.0000001,
				"Description One", 10, 100, 10, 100, 200, 696.0, "swissprot");
		br.setDomainSimilarityScore(0.9);
		getSettings().setDescriptionScoreDomainSimilarityWeight(0.5);

		Double dss = p.getDescriptionScoreCalculator()
				.domainSimilarityScore(br);
		assertNotNull(
				"The description score domain similarity score should have been calculated.",
				dss);
		assertEquals(0.45, dss, 0.001);
	}
}
