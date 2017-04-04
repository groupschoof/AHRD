package ahrd.test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.assertNotNull;
import static ahrd.controller.Settings.getSettings;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.junit.Before;
import org.junit.Test;

import ahrd.model.CompetitorAnnotation;
import ahrd.model.BlastResult;
import ahrd.model.EvaluationScoreCalculator;
import ahrd.model.Fscore;
import ahrd.model.GOterm;
import ahrd.model.Protein;
import ahrd.model.GroundTruthDescription;
import ahrd.model.GOdatabase;

public class EvaluationScoreCalculatorTest {
	private Map<String, GOterm> goDB;
	
	@Before
	public void setup() throws IOException {
		// Init:
		TestUtils.initTestSettings();
	}

	@Test
	public void testTruePositives() {
		Set<String> groundTruthTokens = new HashSet<String>(Arrays.asList(
				"sheep", "wool", "growth", "factor"));
		Set<String> assignedDesc1 = new HashSet<String>(Arrays.asList("goat",
				"wool", "growth"));
		Set<String> assignedDesc2 = new HashSet<String>();
		// Should have 2 TP:
		assertEquals(2.0, EvaluationScoreCalculator.truePositives(
				assignedDesc1, groundTruthTokens), 0.0);
		// Should have zero TP:
		assertEquals(0.0, EvaluationScoreCalculator.truePositives(
				assignedDesc2, groundTruthTokens), 0.0);
	}

	@Test
	public void testTruePositivesRate() {
		Set<String> groundTruthTokens = new HashSet<String>(Arrays.asList(
				"sheep", "wool", "growth", "factor"));
		Set<String> assignedDesc1 = new HashSet<String>(Arrays.asList("goat",
				"wool", "growth"));
		Set<String> assignedDesc2 = new HashSet<String>();
		// Should be 2/4
		assertEquals(0.5, EvaluationScoreCalculator.fBetaScore(assignedDesc1, groundTruthTokens).getRecall(), 0.0);
		// Should be 4/4
		assertEquals(1.0, EvaluationScoreCalculator.fBetaScore(groundTruthTokens, groundTruthTokens).getRecall(), 0.0);
		// Should be 0/4
		assertEquals(0.0, EvaluationScoreCalculator.fBetaScore(assignedDesc2, groundTruthTokens).getRecall(), 0.0);
	}

	@Test
	public void testFalsePositivesRate() {
		Set<String> allBlastTokens = new HashSet<String>(Arrays.asList("sheep",
				"wool", "growth", "factor", "goat", "horn", "tail"));
		Set<String> groundTruthTokens = new HashSet<String>(Arrays.asList(
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
				assignedDesc1, groundTruthTokens, allBlastTokens), 0.0);
		// Should be 0 / 3
		assertEquals(0.0, EvaluationScoreCalculator.falsePositivesRate(
				groundTruthTokens, groundTruthTokens, allBlastTokens), 0.0);
		// Should be 0 / 3
		assertEquals(0.0, EvaluationScoreCalculator.falsePositivesRate(
				assignedDesc2, groundTruthTokens, allBlastTokens), 0.0);
		// Should be 3 / 3
		assertEquals(1.0, EvaluationScoreCalculator.falsePositivesRate(
				assignedDesc3, groundTruthTokens, allBlastTokens), 0.0);
	}

	@Test
	public void testF1Score() {
		Set<String> groundTruthTokens = new HashSet<String>(Arrays.asList(
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
				groundTruthTokens).getScore(), 0.0);
		// 2 * 2/4 * 2/4 / (2/4 + 2/4) = 0.5
		assertEquals(0.5, EvaluationScoreCalculator.fBetaScore(assignedDesc2,
				groundTruthTokens).getScore(), 0.0);
		// 2 * 2/4 * 2/2 / (2/4 + 2/2) = 2/3
		assertEquals((2.0 / 3.0), EvaluationScoreCalculator.fBetaScore(
				assignedDesc3, groundTruthTokens).getScore(), 0.0);
		// The order of tokens shouldn't have any effect:
		assertEquals(1.0, EvaluationScoreCalculator.fBetaScore(assignedDesc4,
				groundTruthTokens).getScore(), 0.0);
		// Zero true-positives should result in a f1-score of zero:
		assertEquals(0.0, EvaluationScoreCalculator.fBetaScore(assignedDesc5,
				groundTruthTokens).getScore(), 0.0);
	}

	@Test
	public void testAssignEvlScrsToCompetitors() throws FileNotFoundException, IOException {
		Protein p = TestUtils.mockProtein();
		// Mock ground truth description:
		GroundTruthDescription rd = new GroundTruthDescription();
		rd.setAccession("AHRDv2_Acc");
		rd.setDescription("AHRD is the best annotator");
		rd.setTokens(new HashSet<String>(Arrays.asList("ahrd", "is", "the",
				"best", "annotator")));
		p.getEvaluationScoreCalculator().setGroundTruthDescription(rd);
		// mock description assigned by AHRD:
		BlastResult ahrdsRes = TestUtils.mockBlastResult(
				"AHRDv2_Acc",
				0.001,
				"AHRD is the best annotator",
				0,
				200,
				10,
				20,
				200,
				30000.0,
				"swissprot",
				new HashSet<String>(Arrays.asList("ahrd", "is", "the", "best",
						"annotator")));
		p.getDescriptionScoreCalculator()
				.setHighestScoringBlastResult(ahrdsRes);
		// mock competitive results
		p.getEvaluationScoreCalculator().addBestUnchangedBlastResult(
				"swissprot",
				TestUtils.mockBlastResult(
						"Sprot One",
						0.001,
						"AHRD is the best",
						0,
						200,
						10,
						20,
						200,
						30000.0,
						"swissprot",
						new HashSet<String>(Arrays.asList("ahrd", "is", "the",
								"best"))));
		p.getEvaluationScoreCalculator().addBestUnchangedBlastResult(
				"trembl",
				TestUtils.mockBlastResult(
						"trEMBL One",
						0.001,
						"AHRD is best eaten alive",
						0,
						200,
						10,
						20,
						200,
						30000.0,
						"trembl",
						new HashSet<String>(Arrays.asList("ahrd", "is", "best",
								"eaten", "alive"))));
		p.getEvaluationScoreCalculator().addBestUnchangedBlastResult(
				"tair",
				TestUtils.mockBlastResult(
						"TAIR One",
						0.001,
						"AHRD is a sheep",
						0,
						200,
						10,
						20,
						200,
						30000.0,
						"tair",
						new HashSet<String>(Arrays.asList("ahrd", "is", "a",
								"sheep"))));
		// mock CompetitorAnnotations:
		if (goDB == null) {
			goDB = new GOdatabase().getMap();
		}
		CompetitorAnnotation firstAnnot = new CompetitorAnnotation("AHRDv2_Acc", "AHRD horn growthase");
		firstAnnot.setGoAnnotations(new HashSet<GOterm>(Arrays.asList(goDB.get("GO:0005524"))));
		firstAnnot.setEvaluationScore(new Fscore(0.1, 0.0, 0.0));
		CompetitorAnnotation secondAnnot = new CompetitorAnnotation("AHRDv2_Acc", "Goat wool growthase");
		secondAnnot.setGoAnnotations(new HashSet<GOterm>(Arrays.asList(goDB.get("GO:0009853"))));
		secondAnnot.setEvaluationScore(new Fscore(0.5, 0.0, 0.0));
		p.getEvaluationScoreCalculator().addCompetitorAnnotation("blast2go", firstAnnot);
		p.getEvaluationScoreCalculator().addCompetitorAnnotation("eggNOGmapper", secondAnnot);
		Map<String, String> blast2goSettings = new HashMap<String, String>();
		blast2goSettings.put("descriptions:", "./mock_file_path.tsv");
		blast2goSettings.put("go_annotations:", "./mock_file_path.goa");
		Map<String, String> eggNOGmapperSettings = new HashMap<String, String>();
		eggNOGmapperSettings.put("descriptions:", "./mock_file_path.tsv");
		eggNOGmapperSettings.put("go_annotations:", "./mock_file_path.goa");
		Map<String, Map<String, String>> competitorSettings = new HashMap<String, Map<String, String>>();
		competitorSettings.put("blast2go", blast2goSettings);
		competitorSettings.put("eggNOGmapper", eggNOGmapperSettings);
		getSettings().setCompetitorSettings(competitorSettings);
		// TEST:
		p.getEvaluationScoreCalculator().assignEvaluationScores();
		// Test Assignment of scores to the three best Blast-Hits (sprot, tair,
		// trembl):
		assertEquals(0.888888888888889, p.getEvaluationScoreCalculator()
				.getBestUnchangedBlastResults().get("swissprot")
				.getEvaluationScore().getScore(), 0.0);
		assertEquals(0.6, p.getEvaluationScoreCalculator()
				.getBestUnchangedBlastResults().get("trembl").getEvaluationScore().getScore(),
				0.0);
		assertEquals(0.4444444444444445, p.getEvaluationScoreCalculator()
				.getBestUnchangedBlastResults().get("tair").getEvaluationScore().getScore(),
				0.0);
		// Test Assignment of scores to the two CompetitorAnnotations:
		assertEquals(0.25, firstAnnot.getEvaluationScore().getScore(), 0.0);
		assertEquals(0.0, secondAnnot.getEvaluationScore().getScore(), 0.0);
		// evaluation-scores should be fine (tested above)
		// Just test the resulting evalScoreMinBestCompScore:
		assertEquals(0.111111111111111, p.getEvaluationScoreCalculator()
				.getEvalScoreMinBestCompScore(), 0.000001);
	}

	@Test
	public void testFindHighestPossibleEvaluationScore() {
		// Mock Protein:
		Protein p = TestUtils.mockProtein();
		// Mock ground truth description:
		GroundTruthDescription rd = new GroundTruthDescription();
		rd.setAccession("AHRDv2_Acc");
		rd.setDescription("AHRD is the best annotator");
		rd.setTokens(new HashSet<String>(Arrays.asList("ahrd", "is", "the",
				"best", "annotator")));
		p.getEvaluationScoreCalculator().setGroundTruthDescription(rd);
		// Mock first BlastResults:
		List<BlastResult> swissprotRes = new ArrayList<BlastResult>();
		swissprotRes.add(TestUtils.mockBlastResult(
				"AHRDv2_Acc",
				0.001,
				"AHRD is the best annotator",
				0,
				200,
				10,
				20,
				200,
				30000.0,
				"swissprot",
				new HashSet<String>(Arrays.asList("ahrd", "is", "the", "best",
						"annotator"))));
		// ... and the second BlastResult:
		swissprotRes.add(TestUtils
				.mockBlastResult(
						"Sprot One",
						0.001,
						"AHRD is the best",
						0,
						200,
						10,
						20,
						200,
						30000.0,
						"swissprot",
						new HashSet<String>(Arrays.asList("ahrd", "is", "the",
								"best"))));
		p.getBlastResults().put("swissprot", swissprotRes);
		// test:
		p.getEvaluationScoreCalculator().findBlastResultWithHighestPossibleDescriptionScore();
		assertNotNull("Highest possible evaluation score should be not null", p
				.getEvaluationScoreCalculator()
				.getHighestPossibleDescriptionScore());
		assertEquals("Highest possible evaluation score should be 1.0", p
				.getEvaluationScoreCalculator()
				.getHighestPossibleDescriptionScore().getScore(), 1.0, 0.0);
	}

	@Test
	public void testAddCompetitorAnnotion() throws FileNotFoundException, IOException {
		Protein p = TestUtils.mockProtein();
		if (goDB == null) {
			goDB = new GOdatabase().getMap();
		}
		CompetitorAnnotation firstAnnot = new CompetitorAnnotation("P04637", "Cellular tumor antigen p53");
		firstAnnot.setGoAnnotations(new HashSet<GOterm>(Arrays.asList(goDB.get("GO:0005524"))));
		firstAnnot.setEvaluationScore(new Fscore(0.1, 0.0, 0.0));
		CompetitorAnnotation firstAnnotClone = new CompetitorAnnotation("P04637", "Cellular tumor antigen p53");
		firstAnnotClone.setGoAnnotations(new HashSet<GOterm>(Arrays.asList(goDB.get("GO:0005524"))));
		firstAnnotClone.setEvaluationScore(new Fscore(0.1, 0.0, 0.0));
		CompetitorAnnotation secondAnnot = new CompetitorAnnotation("P0C512", "Ribulose bisphosphate carboxylase large chain");
		secondAnnot.setGoAnnotations(new HashSet<GOterm>(Arrays.asList(goDB.get("GO:0009853"))));
		// test
		p.getEvaluationScoreCalculator().addCompetitorAnnotation("blast2go", firstAnnot);
		assertEquals(1, p.getEvaluationScoreCalculator().getCompetitorAnnotations().size());
		assertEquals("Cellular tumor antigen p53", p.getEvaluationScoreCalculator().getCompetitorAnnotations().get("blast2go").getDescription());
		p.getEvaluationScoreCalculator().addCompetitorAnnotation("blast2go", firstAnnotClone);
		assertEquals(1, p.getEvaluationScoreCalculator().getCompetitorAnnotations().size());
		assertEquals("Cellular tumor antigen p53", p.getEvaluationScoreCalculator().getCompetitorAnnotations().get("blast2go").getDescription());
		p.getEvaluationScoreCalculator().addCompetitorAnnotation("eggNOG-mapper", secondAnnot);
		assertEquals(2, p.getEvaluationScoreCalculator().getCompetitorAnnotations().size());
		for (CompetitorAnnotation compAnnot : p.getEvaluationScoreCalculator().getCompetitorAnnotations().values()) {
			assertTrue("Cellular tumor antigen p53".equals(compAnnot.getDescription())	|| "Ribulose bisphosphate carboxylase large chain".equals(compAnnot.getDescription()));
		}
	}

	@Test
	public void testCalcSimpleGoAnnotationScore() {
		getSettings().setPathToGeneOntologyReference("swissprot","./test/resources/reference_gene_ontology_annotations_uniprotKB_GOA.txt");
		getSettings().setPathToGroundTruthGoAnnotations("./test/resources/sprot_GO_ground_truth.goa");
		getSettings().setCalculateSimpleGoF1Scores(true);
		Protein p = TestUtils.mockProtein();
		p.setEvaluationScoreCalculator(new EvaluationScoreCalculator(p));
		GOterm bpRoot = new GOterm("GO:0008150", "biological_process", "biological_process");
		GOterm bpCellularProcess = new GOterm("GO:0009987", "cellular process", "biological_process");
		GOterm bpSingleOrganismProcess = new GOterm("GO:0044699", "single-organism process", "biological_process");
		// |ref| == 0 && |pred| == 0 -> f1 == 1
		p.getEvaluationScoreCalculator().assignEvaluationScores();
		assertTrue(p.getEvaluationScoreCalculator().getSimpleGoAnnotationScore().getScore() == 1.0);
		// |ref| == 0 && |pred| > 0 -> f1 == 0
		p.getEvaluationScoreCalculator().setGroundTruthGoAnnoatations(new HashSet<GOterm>(Arrays.asList(bpRoot)));
		p.getEvaluationScoreCalculator().assignEvaluationScores();
		assertTrue(p.getEvaluationScoreCalculator().getSimpleGoAnnotationScore().getScore() == 0.0);
		// |ref| > 0 && |pred| == 0 -> f1 == 0
		p.getEvaluationScoreCalculator().setGroundTruthGoAnnoatations(new HashSet<GOterm>());
		p.setGoResultsTerms(new HashSet<GOterm>(Arrays.asList(bpRoot)));
		p.getEvaluationScoreCalculator().assignEvaluationScores();
		assertTrue(p.getEvaluationScoreCalculator().getSimpleGoAnnotationScore().getScore() == 0.0);
		// ref == pred -> f1 == 1
		p.getEvaluationScoreCalculator().setGroundTruthGoAnnoatations(new HashSet<GOterm>(Arrays.asList(bpRoot, bpCellularProcess, bpSingleOrganismProcess)));
		p.setGoResultsTerms(new HashSet<GOterm>(Arrays.asList(bpRoot, bpCellularProcess, bpSingleOrganismProcess)));
		p.getEvaluationScoreCalculator().assignEvaluationScores();
		assertTrue(p.getEvaluationScoreCalculator().getSimpleGoAnnotationScore().getScore() == 1.0);
		// |ref| > 0 && |pred| > 0 && ref != pred -> f1 == 0
		p.getEvaluationScoreCalculator().setGroundTruthGoAnnoatations(new HashSet<GOterm>(Arrays.asList(bpCellularProcess)));
		p.setGoResultsTerms(new HashSet<GOterm>(Arrays.asList(bpSingleOrganismProcess)));
		p.getEvaluationScoreCalculator().assignEvaluationScores();
		assertTrue(p.getEvaluationScoreCalculator().getSimpleGoAnnotationScore().getScore() == 0.0);
		// |ref| == 2 && |pred| == 2 && |ref^pred| == 1 -> f1 == 2*0.5*0.5/(0.5+0.5) == 0.5
		p.getEvaluationScoreCalculator().setGroundTruthGoAnnoatations(new HashSet<GOterm>(Arrays.asList(bpRoot, bpCellularProcess)));
		p.setGoResultsTerms(new HashSet<GOterm>(Arrays.asList(bpRoot, bpSingleOrganismProcess)));
		p.getEvaluationScoreCalculator().assignEvaluationScores();
		assertTrue(p.getEvaluationScoreCalculator().getSimpleGoAnnotationScore().getScore() == 0.5);
		// |ref| == 1 && |pred| == 2 && |ref^pred| == 1 -> f1 == 2*1*0.5/(1+0.5) == 0.667
		p.getEvaluationScoreCalculator().setGroundTruthGoAnnoatations(new HashSet<GOterm>(Arrays.asList(bpCellularProcess)));
		p.setGoResultsTerms(new HashSet<GOterm>(Arrays.asList(bpCellularProcess, bpSingleOrganismProcess)));
		p.getEvaluationScoreCalculator().assignEvaluationScores();
		assertTrue(p.getEvaluationScoreCalculator().getSimpleGoAnnotationScore().getScore() == (double)2/3);
	}
	
	@Test
	public void testCalcAncestryGoAnnotationScore() {
		getSettings().setPathToGeneOntologyReference("swissprot","./test/resources/reference_gene_ontology_annotations_uniprotKB_GOA.txt");
		getSettings().setPathToGroundTruthGoAnnotations("./test/resources/sprot_GO_ground_truth.goa");
		getSettings().setCalculateAncestryGoF1Scores(true);
		Protein p = TestUtils.mockProtein();
		p.setEvaluationScoreCalculator(new EvaluationScoreCalculator(p));
		GOterm bpRoot = new GOterm("GO:0008150", "biological_process", "biological_process");
		bpRoot.addTermToAncestry(bpRoot);
		GOterm bpCellularProcess = new GOterm("GO:0009987", "cellular process", "biological_process");
		bpCellularProcess.addTermToAncestry(bpRoot);
		bpCellularProcess.addTermToAncestry(bpCellularProcess);
		GOterm bpSingleOrganismProcess = new GOterm("GO:0044699", "single-organism process", "biological_process");
		bpSingleOrganismProcess.addTermToAncestry(bpRoot);
		bpSingleOrganismProcess.addTermToAncestry(bpSingleOrganismProcess);
		// |ref| == 0 && |pred| == 0 -> f1 == 1
		p.getEvaluationScoreCalculator().assignEvaluationScores();
		assertTrue(p.getEvaluationScoreCalculator().getAncestryGoAnnotationScore().getScore() == 1.0);
		// |ref| == 0 && |pred| > 0 -> f1 == 0
		p.getEvaluationScoreCalculator().setGroundTruthGoAnnoatations(new HashSet<GOterm>(Arrays.asList(bpRoot)));
		p.getEvaluationScoreCalculator().assignEvaluationScores();
		assertTrue(p.getEvaluationScoreCalculator().getAncestryGoAnnotationScore().getScore() == 0.0);
		// |ref| > 0 && |pred| == 0 -> f1 == 0
		p.getEvaluationScoreCalculator().setGroundTruthGoAnnoatations(new HashSet<GOterm>());
		p.setGoResultsTerms(new HashSet<GOterm>(Arrays.asList(bpRoot)));
		p.getEvaluationScoreCalculator().assignEvaluationScores();
		assertTrue(p.getEvaluationScoreCalculator().getAncestryGoAnnotationScore().getScore() == 0.0);
		// ref == pred -> f1 == 1
		p.getEvaluationScoreCalculator().setGroundTruthGoAnnoatations(new HashSet<GOterm>(Arrays.asList(bpRoot, bpCellularProcess, bpSingleOrganismProcess)));
		p.setGoResultsTerms(new HashSet<GOterm>(Arrays.asList(bpRoot, bpCellularProcess, bpSingleOrganismProcess)));
		p.getEvaluationScoreCalculator().assignEvaluationScores();
		assertTrue(p.getEvaluationScoreCalculator().getAncestryGoAnnotationScore().getScore() == 1.0);
		// |ref| == 1 && |pred| == 1 && ref != pred && |ancestry(ref)| == |ancestry(pred)| == 2 && |ancestry(ref)^ancestry(pred)| == 1 -> f1 == 2*0.5*0.5/(0.5+0.5) == 0.5
		p.getEvaluationScoreCalculator().setGroundTruthGoAnnoatations(new HashSet<GOterm>(Arrays.asList(bpCellularProcess)));
		p.setGoResultsTerms(new HashSet<GOterm>(Arrays.asList(bpSingleOrganismProcess)));
		p.getEvaluationScoreCalculator().assignEvaluationScores();
		assertTrue(p.getEvaluationScoreCalculator().getAncestryGoAnnotationScore().getScore() == 0.5);
		// ref == ancestry(ref) && pred == ancestry(pred) && |ref| == 2 && |pred| == 2 && ref != pred && |ancestry(ref)| == |ancestry(pred)| == 2 && |ancestry(ref)^ancestry(pred)| == 1 -> f1 == 2*0.5*0.5/(0.5+0.5) == 0.5
		p.getEvaluationScoreCalculator().setGroundTruthGoAnnoatations(new HashSet<GOterm>(Arrays.asList(bpRoot, bpCellularProcess)));
		p.setGoResultsTerms(new HashSet<GOterm>(Arrays.asList(bpRoot, bpSingleOrganismProcess)));
		p.getEvaluationScoreCalculator().assignEvaluationScores();
		assertTrue(p.getEvaluationScoreCalculator().getAncestryGoAnnotationScore().getScore() == 0.5);
	}
	@Test
	public void testCalcSemSimGoAnnotationScore() {
		getSettings().setPathToGeneOntologyReference("swissprot","./test/resources/reference_gene_ontology_annotations_uniprotKB_GOA.txt");
		getSettings().setPathToGroundTruthGoAnnotations("./test/resources/sprot_GO_ground_truth.goa");
		getSettings().setCalculateSemSimGoF1Scores(true);
		Protein p = TestUtils.mockProtein();
		p.setEvaluationScoreCalculator(new EvaluationScoreCalculator(p));
		GOterm mfRoot = new GOterm("GO:0003674", "molecular_function", "molecular_function");
		mfRoot.addTermToAncestry(mfRoot);
		mfRoot.setInformationContent(0.0);
		GOterm bpRoot = new GOterm("GO:0008150", "biological_process", "biological_process");
		bpRoot.addTermToAncestry(bpRoot);
		bpRoot.setInformationContent(0.0);
		GOterm bpCellularProcess = new GOterm("GO:0009987", "cellular process", "biological_process");
		bpCellularProcess.addTermToAncestry(bpRoot);
		bpCellularProcess.addTermToAncestry(bpCellularProcess);
		bpCellularProcess.setInformationContent(0.423452947461567);
		GOterm bpSingleOrganismProcess = new GOterm("GO:0044699", "single-organism process", "biological_process");
		bpSingleOrganismProcess.addTermToAncestry(bpRoot);
		bpSingleOrganismProcess.addTermToAncestry(bpSingleOrganismProcess);
		bpSingleOrganismProcess.setInformationContent(0.734371256710225);
		// |ref| == 0 && |pred| == 0 -> f1 == 1
		p.getEvaluationScoreCalculator().assignEvaluationScores();
		assertTrue(p.getEvaluationScoreCalculator().getSemSimGoAnnotationScore().getScore() == 1.0);
		// |ref| == 0 && |pred| > 0 -> f1 == 0
		p.getEvaluationScoreCalculator().setGroundTruthGoAnnoatations(new HashSet<GOterm>(Arrays.asList(bpRoot)));
		p.getEvaluationScoreCalculator().assignEvaluationScores();
		assertTrue(p.getEvaluationScoreCalculator().getSemSimGoAnnotationScore().getScore() == 0.0);
		// |ref| > 0 && |pred| == 0 -> f1 == 0
		p.getEvaluationScoreCalculator().setGroundTruthGoAnnoatations(new HashSet<GOterm>());
		p.setGoResultsTerms(new HashSet<GOterm>(Arrays.asList(bpRoot)));
		p.getEvaluationScoreCalculator().assignEvaluationScores();
		assertTrue(p.getEvaluationScoreCalculator().getSemSimGoAnnotationScore().getScore() == 0.0);
		// ref == pred -> f1 == 1
		p.getEvaluationScoreCalculator().setGroundTruthGoAnnoatations(new HashSet<GOterm>(Arrays.asList(bpCellularProcess, bpSingleOrganismProcess)));
		p.setGoResultsTerms(new HashSet<GOterm>(Arrays.asList(bpCellularProcess, bpSingleOrganismProcess)));
		p.getEvaluationScoreCalculator().assignEvaluationScores();
		assertTrue(p.getEvaluationScoreCalculator().getSemSimGoAnnotationScore().getScore() == 1.0);
		// |ref| == 1 && |pred| == 1 && ref != pred && commonAncestor(ref, pred) == rootTerm -> f1 == 0
		p.getEvaluationScoreCalculator().setGroundTruthGoAnnoatations(new HashSet<GOterm>(Arrays.asList(bpCellularProcess)));
		p.setGoResultsTerms(new HashSet<GOterm>(Arrays.asList(bpSingleOrganismProcess)));
		p.getEvaluationScoreCalculator().assignEvaluationScores();
		assertTrue(p.getEvaluationScoreCalculator().getSemSimGoAnnotationScore().getScore() == 0.0);
		// ref == root && pred != root -> f1 == 0
		p.getEvaluationScoreCalculator().setGroundTruthGoAnnoatations(new HashSet<GOterm>(Arrays.asList(bpRoot)));
		p.setGoResultsTerms(new HashSet<GOterm>(Arrays.asList(bpSingleOrganismProcess)));
		p.getEvaluationScoreCalculator().assignEvaluationScores();
		assertTrue(p.getEvaluationScoreCalculator().getSemSimGoAnnotationScore().getScore() == 0.0);
		// ref == pred == root -> f1 == 1
		p.getEvaluationScoreCalculator().setGroundTruthGoAnnoatations(new HashSet<GOterm>(Arrays.asList(bpRoot)));
		p.setGoResultsTerms(new HashSet<GOterm>(Arrays.asList(bpRoot)));
		p.getEvaluationScoreCalculator().assignEvaluationScores();
		assertTrue(p.getEvaluationScoreCalculator().getSemSimGoAnnotationScore().getScore() == 1.0);
		// ref == root && pred == root && ref != pred -> f1 == 1
		p.getEvaluationScoreCalculator().setGroundTruthGoAnnoatations(new HashSet<GOterm>(Arrays.asList(mfRoot)));
		p.setGoResultsTerms(new HashSet<GOterm>(Arrays.asList(bpRoot)));
		p.getEvaluationScoreCalculator().assignEvaluationScores();
		assertTrue(p.getEvaluationScoreCalculator().getSemSimGoAnnotationScore().getScore() == 1.0);
	}
	
}
