package ahrd.test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.junit.Before;
import org.junit.Test;

import ahrd.model.Blast2GoAnnot;
import ahrd.model.BlastResult;
import ahrd.model.EvaluationScoreCalculator;
import ahrd.model.Protein;
import ahrd.model.ReferenceDescription;

public class EvaluationScoreCalculatorTest {

	@Before
	public void setup() throws IOException {
		// Init:
		TestUtils.initTestSettings();
	}

	@Test
	public void testTruePositives() {
		Set<String> referenceTokens = new HashSet<String>(Arrays.asList(
				"sheep", "wool", "growth", "factor"));
		Set<String> assignedDesc1 = new HashSet<String>(Arrays.asList("goat",
				"wool", "growth"));
		Set<String> assignedDesc2 = new HashSet<String>();
		// Should have 2 TP:
		assertEquals(2.0, EvaluationScoreCalculator.truePositives(
				assignedDesc1, referenceTokens), 0.0);
		// Should have zero TP:
		assertEquals(0.0, EvaluationScoreCalculator.truePositives(
				assignedDesc2, referenceTokens), 0.0);
	}

	@Test
	public void testTruePositivesRate() {
		Set<String> referenceTokens = new HashSet<String>(Arrays.asList(
				"sheep", "wool", "growth", "factor"));
		Set<String> assignedDesc1 = new HashSet<String>(Arrays.asList("goat",
				"wool", "growth"));
		Set<String> assignedDesc2 = new HashSet<String>();
		// Should be 2/4
		assertEquals(0.5, EvaluationScoreCalculator.truePositivesRate(
				assignedDesc1, referenceTokens), 0.0);
		// Should be 4/4
		assertEquals(1.0, EvaluationScoreCalculator.truePositivesRate(
				referenceTokens, referenceTokens), 0.0);
		// Should be 0/4
		assertEquals(0.0, EvaluationScoreCalculator.truePositivesRate(
				assignedDesc2, referenceTokens), 0.0);
	}

	@Test
	public void testFalsePositivesRate() {
		Set<String> allBlastTokens = new HashSet<String>(Arrays.asList("sheep",
				"wool", "growth", "factor", "goat", "horn", "tail"));
		Set<String> referenceTokens = new HashSet<String>(Arrays.asList(
				"sheep", "wool", "growth", "factor"));
		Set<String> assignedDesc1 = new HashSet<String>(Arrays.asList("goat",
				"wool", "growth"));
		Set<String> assignedDesc2 = new HashSet<String>();
		Set<String> assignedDesc3 = new HashSet<String>(Arrays.asList("goat",
				"wool", "growth", "horn", "tail"));

		// TEST:
		// The Set of all negatives has cardinality of 3.

		// Should be 1 / 3
		assertEquals((1.0 / 3.0), EvaluationScoreCalculator.falsePositivesRate(
				assignedDesc1, referenceTokens, allBlastTokens), 0.0);
		// Should be 0 / 3
		assertEquals(0.0, EvaluationScoreCalculator.falsePositivesRate(
				referenceTokens, referenceTokens, allBlastTokens), 0.0);
		// Should be 0 / 3
		assertEquals(0.0, EvaluationScoreCalculator.falsePositivesRate(
				assignedDesc2, referenceTokens, allBlastTokens), 0.0);
		// Should be 3 / 3
		assertEquals(1.0, EvaluationScoreCalculator.falsePositivesRate(
				assignedDesc3, referenceTokens, allBlastTokens), 0.0);
	}

	@Test
	public void testF1Score() {
		Set<String> referenceTokens = new HashSet<String>(Arrays.asList(
				"sheep", "wool", "growth", "factor"));
		Set<String> assignedDesc1 = new HashSet<String>(Arrays.asList("sheep",
				"wool", "growth", "factor"));
		Set<String> assignedDesc2 = new HashSet<String>(Arrays.asList("goat",
				"horn", "growth", "factor"));
		Set<String> assignedDesc3 = new HashSet<String>(Arrays.asList("sheep",
				"wool"));
		Set<String> assignedDesc4 = new HashSet<String>(Arrays.asList("growth",
				"sheep", "factor", "wool"));
		Set<String> assignedDesc5 = new HashSet<String>(Arrays.asList("true",
				"goats", "do", "not", "share", "tokens"));

		// F1-score is the harmonic mean of precision and recall, where
		// precision := #true-positives / #assigned-tokens and
		// recall := #true-positives / #reference-tokens.

		// 2 * 4/4 * 4/4 / (4/4 + 4/4) = 1
		assertEquals(1.0, EvaluationScoreCalculator.fBetaScore(assignedDesc1,
				referenceTokens), 0.0);
		// 2 * 2/4 * 2/4 / (2/4 + 2/4) = 0.5
		assertEquals(0.5, EvaluationScoreCalculator.fBetaScore(assignedDesc2,
				referenceTokens), 0.0);
		// 2 * 2/4 * 2/2 / (2/4 + 2/2) = 2/3
		assertEquals((2.0 / 3.0), EvaluationScoreCalculator.fBetaScore(
				assignedDesc3, referenceTokens), 0.0);
		// The order of tokens shouldn't have any effect:
		assertEquals(1.0, EvaluationScoreCalculator.fBetaScore(assignedDesc4,
				referenceTokens), 0.0);
		// Zero true-positives should result in a f1-score of zero:
		assertEquals(0.0, EvaluationScoreCalculator.fBetaScore(assignedDesc5,
				referenceTokens), 0.0);
	}

	@Test
	public void testAssignEvlScrsToCompetitors() {
		Protein p = TestUtils.mockProtein();
		// Mock reference description:
		ReferenceDescription rd = new ReferenceDescription();
		rd.setAccession("AHRDv2_Acc");
		rd.setDescription("AHRD is the best annotator");
		rd.setTokens(new HashSet<String>(Arrays.asList("ahrd", "is", "the",
				"best", "annotator")));
		p.getEvaluationScoreCalculator().setReferenceDescription(rd);
		// mock description assigned by AHRD:
		BlastResult ahrdsRes = TestUtils.mockBlastResult("AHRDv2_Acc", 0.001,
				"AHRD is the best annotator", 0, 200, 30000.0, "swissprot",
				new HashSet<String>(Arrays.asList("ahrd", "is", "the", "best",
						"annotator")));
		p.getDescriptionScoreCalculator()
				.setHighestScoringBlastResult(ahrdsRes);
		// mock competitive results
		p.getEvaluationScoreCalculator().addUnchangedBlastResult(
				"swissprot",
				TestUtils.mockBlastResult("Sprot One", 0.001,
						"AHRD is the best", 0, 200, 30000.0, "swissprot",
						new HashSet<String>(Arrays.asList("ahrd", "is", "the",
								"best"))));
		p.getEvaluationScoreCalculator().addUnchangedBlastResult(
				"trembl",
				TestUtils.mockBlastResult("trEMBL One", 0.001,
						"AHRD is best eaten alive", 0, 200, 30000.0, "trembl",
						new HashSet<String>(Arrays.asList("ahrd", "is", "best",
								"eaten", "alive"))));
		p.getEvaluationScoreCalculator().addUnchangedBlastResult(
				"tair",
				TestUtils.mockBlastResult("TAIR One", 0.001, "AHRD is a sheep",
						0, 200, 30000.0, "tair", new HashSet<String>(Arrays
								.asList("ahrd", "is", "a", "sheep"))));
		// mock Blast2GoAnnots:
		String b2gResEntry = "accession\tGO:123456\tAHRD horn growthase";
		String b2gResEntryTwo = "accession\tGO:654321\tGoat wool growthase";
		Blast2GoAnnot b2gaRef = Blast2GoAnnot.fromBlast2GoEntry(b2gResEntry);
		b2gaRef.setEvaluationScore(0.1);
		Blast2GoAnnot b2gaTwo = Blast2GoAnnot.fromBlast2GoEntry(b2gResEntryTwo);
		b2gaTwo.setEvaluationScore(0.5);
		p.getEvaluationScoreCalculator().addBlast2GoAnnot(b2gaRef);
		p.getEvaluationScoreCalculator().addBlast2GoAnnot(b2gaTwo);
		// TEST:
		p.getEvaluationScoreCalculator().assignEvlScrsToCompetitors();
		// Test Assignment of scores to the three best Blast-Hits (sprot, tair,
		// trembl):
		assertEquals(0.888888888888889, p.getEvaluationScoreCalculator()
				.getUnchangedBlastResults().get("swissprot")
				.getEvaluationScore(), 0.0);
		assertEquals(0.6, p.getEvaluationScoreCalculator()
				.getUnchangedBlastResults().get("trembl").getEvaluationScore(),
				0.0);
		assertEquals(0.4444444444444445, p.getEvaluationScoreCalculator()
				.getUnchangedBlastResults().get("tair").getEvaluationScore(),
				0.0);
		// Test Assignment of scores to the two Blat2GoAnnots:
		assertEquals(0.25, b2gaRef.getEvaluationScore(), 0.0);
		assertEquals(0.0, b2gaTwo.getEvaluationScore(), 0.0);
		// evaluation-scores should be fine (tested above)
		// Just test the resulting evalScoreMinBestCompScore:
		assertEquals(0.111111111111111, p.getEvaluationScoreCalculator()
				.getEvalScoreMinBestCompScore(), 0.000001);

		// test, if TPR and FPR have been set:
		assertTrue("FPR should have been set, but is NULL.", p
				.getEvaluationScoreCalculator().getFalsePositivesRate() != null);
		assertTrue("TPR should have been set, but is NULL.", p
				.getEvaluationScoreCalculator().getTruePositivesRate() != null);
	}

	@Test
	public void testAddBlast2GoAnnot() {
		Protein p = TestUtils.mockProtein();
		String b2gResEntry = "accession\tGO:123456\tSheep horn growthase";
		String b2gResEntryTwo = "accession\tGO:654321\tGoat wool grothase";
		Blast2GoAnnot b2gaRef = Blast2GoAnnot.fromBlast2GoEntry(b2gResEntry);
		Blast2GoAnnot b2gaClone = Blast2GoAnnot.fromBlast2GoEntry(b2gResEntry);
		Blast2GoAnnot b2gaTwo = Blast2GoAnnot.fromBlast2GoEntry(b2gResEntryTwo);
		// test
		p.getEvaluationScoreCalculator().addBlast2GoAnnot(b2gaRef);
		assertEquals(1, p.getEvaluationScoreCalculator().getBlast2GoAnnots()
				.size());
		assertEquals("Sheep horn growthase", p.getEvaluationScoreCalculator()
				.getBlast2GoAnnots().toArray(new Blast2GoAnnot[] {})[0]
				.getDescription());
		p.getEvaluationScoreCalculator().addBlast2GoAnnot(b2gaClone);
		assertEquals(1, p.getEvaluationScoreCalculator().getBlast2GoAnnots()
				.size());
		assertEquals("Sheep horn growthase", p.getEvaluationScoreCalculator()
				.getBlast2GoAnnots().toArray(new Blast2GoAnnot[] {})[0]
				.getDescription());
		p.getEvaluationScoreCalculator().addBlast2GoAnnot(b2gaTwo);
		assertEquals(2, p.getEvaluationScoreCalculator().getBlast2GoAnnots()
				.size());
		for (Blast2GoAnnot b2ga : p.getEvaluationScoreCalculator()
				.getBlast2GoAnnots().toArray(new Blast2GoAnnot[] {})) {
			assertTrue("Sheep horn growthase".equals(b2ga.getDescription())
					|| "Goat wool grothase".equals(b2ga.getDescription()));
		}
	}

	@Test
	public void testSortBlast2GoAnnotsByEvalScore() {
		Protein p = TestUtils.mockProtein();
		String b2gResEntry = "accession\tGO:123456\tSheep horn growthase";
		String b2gResEntryTwo = "accession\tGO:654321\tGoat wool growthase";
		String b2gResEntryThree = "accession\tGO:654321\tOveja chiquitisima bala mucho";
		Blast2GoAnnot b2gaRef = Blast2GoAnnot.fromBlast2GoEntry(b2gResEntry);
		b2gaRef.setEvaluationScore(0.1);
		Blast2GoAnnot b2gaTwo = Blast2GoAnnot.fromBlast2GoEntry(b2gResEntryTwo);
		b2gaTwo.setEvaluationScore(0.5);
		Blast2GoAnnot b2gaThree = Blast2GoAnnot
				.fromBlast2GoEntry(b2gResEntryThree);
		b2gaThree.setEvaluationScore(0.9);
		Set<Blast2GoAnnot> b2gas = new HashSet<Blast2GoAnnot>();
		b2gas.add(b2gaRef);
		b2gas.add(b2gaTwo);
		b2gas.add(b2gaThree);
		p.getEvaluationScoreCalculator().setBlast2GoAnnots(b2gas);
		// test
		List<Blast2GoAnnot> rankedB2gas = p.getEvaluationScoreCalculator()
				.sortBlast2GoAnnotsByEvalScore();
		assertEquals(b2gaRef, rankedB2gas.get(0));
		assertEquals(b2gaTwo, rankedB2gas.get(1));
		assertEquals(b2gaThree, rankedB2gas.get(2));
	}
	
}
