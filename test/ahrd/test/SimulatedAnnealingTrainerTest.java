package ahrd.test;

import static ahrd.controller.Settings.getSettings;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;
import java.util.logging.Level;

import org.junit.Before;
import org.junit.Test;
import org.semanticweb.owlapi.model.OWLOntologyCreationException;

import ahrd.controller.DescriptionParameters;
import ahrd.controller.Parameters;
import ahrd.controller.SimulatedAnnealingTrainer;
import ahrd.exception.MissingAccessionException;
import ahrd.exception.MissingProteinException;
import ahrd.model.Fscore;
import ahrd.model.Protein;
import ahrd.view.SimulatedAnnealingTrainerOutputWriter;

public class SimulatedAnnealingTrainerTest {

	private SimulatedAnnealingTrainer trainer;

	@Before
	public void setUp() throws IOException, MissingAccessionException, MissingProteinException {
		trainer = new SimulatedAnnealingTrainer("./test/resources/trainer_input.yml");
		getSettings().setLoggingLevel(Level.WARNING);
		trainer.setup();
		trainer.setupGroundTruthDescriptions();
		// Sets up description and GOterm annotations of competitors in the EvaluationScoreCalculators of each Protein
		trainer.setupCompetitors();
	}

	@Test
	public void testAvgEvaluationScore() {
		assertTrue("SimulatedAnnealingTrainer should initialize Settings to Evaluation-Mode.",
				getSettings().isInEvaluationMode());
		Protein p1 = new Protein("protein_one", 200);
		p1.getEvaluationScoreCalculator().setEvalutionScore(new Fscore(1.0, 0.0, 0.6));
		p1.getEvaluationScoreCalculator().setEvalScoreMinBestCompScore(1.0);
		Protein p2 = new Protein("protein_two", 210);
		p2.getEvaluationScoreCalculator().setEvalutionScore(new Fscore(0.8, 0.0, 0.7));
		p2.getEvaluationScoreCalculator().setEvalScoreMinBestCompScore(0.8);
		Protein p3 = new Protein("protein_three", 220);
		p3.getEvaluationScoreCalculator().setEvalutionScore(new Fscore(0.3, 0.0, 0.5));
		p3.getEvaluationScoreCalculator().setEvalScoreMinBestCompScore(0.3);
		this.trainer.setProteins(new HashMap<>());
		this.trainer.getProteins().put(p1.getAccession(), p1);
		this.trainer.getProteins().put(p2.getAccession(), p2);
		this.trainer.getProteins().put(p3.getAccession(), p3);
		// test
		trainer.calcAveragesOfDescriptionScores();
		// (1.0 + 0.8 + 0.3) / 3 = 0.7
		assertEquals(0.7, getSettings().getAvgDescriptionEvaluationScore(), 0.000000000001);
		// (0.6 + 0.7 + 0.5) / 3 = 0.6
		assertEquals(0.6, getSettings().getAvgDescriptionRecall(),
				0.000000000001);
		// test robustness for zero scores:
		for (Protein p : trainer.getProteins().values()) {
			p.getEvaluationScoreCalculator().setEvalScoreMinBestCompScore(0.0);
		}
		trainer.calcAveragesOfDescriptionScores();
		assertEquals(getSettings().getAvgDescriptionEvaluationScore(), 0.7, 0.00000000001);
	}

	@Test
	public void testRememberSimulatedAnnealingPath() throws OWLOntologyCreationException, IOException, MissingAccessionException {
		// Just do two cycles:
		getSettings().setTemperature(2);
		this.trainer.setOutputWriter(new SimulatedAnnealingTrainerOutputWriter(getSettings().getPathToDescriptionTrainingPathLog()));
		Parameters p = getSettings().getDescriptionParameters().clone();
		this.trainer.getOutputWriter().generateHeader(false, p);
		this.trainer.train(p);
		// test:
		assertEquals(
				"Training should have visited two distinct parameter-sets.", 2,
				this.trainer.getTestedParameters().size());
		Parameters clone = this.trainer.getAlreadyTestedParameters(p);
		assertNotNull(clone);
		// Scores should be equal:
		assertNotNull("Avg EvaluationScore should be remembered.",
				clone.getAvgEvaluationScore());
		assertNotNull("Avg Precision (PPV) should be remembered.",
				clone.getAvgPrecision());
		assertNotNull("Avg Recall (TPR) should be remembered.",
				clone.getAvgRecall());
	}

	@Test
	public void testAcceptanceProbability() {
		getSettings().setAvgDescriptionEvaluationScore(0.5);
		getSettings().setOptimizationAcceptanceProbabilityScalingFactor(Double.valueOf(200000000));
		getSettings().setTemperature(1000);
		Parameters p = getSettings().getDescriptionParameters().clone();
		// test first iteration, when accepted Settings are null:
		assertEquals(1.0, trainer.acceptanceProbability(p), 0.0);
		// test current Settings better than accepted:
		trainer.setAcceptedParameters(p);
		getSettings().setAvgDescriptionEvaluationScore(1.0);
		p = getSettings().getDescriptionParameters().clone();
		assertEquals(1.0, trainer.acceptanceProbability(p), 0.0);
		// test current Settings worse than accepted ones:
		trainer.setAcceptedParameters(p);
		getSettings().setAvgDescriptionEvaluationScore(0.9999741);
		p = getSettings().getDescriptionParameters().clone();
		// 0.9999741 - 1.0 = -2.59 * 10^-5
		assertEquals(-0.0000259,
				trainer.diffEvalScoreToCurrentlyAcceptedParams(p), 0.00000001);
		// exp((-0.0000259*200,000,000)/1000) = 0.005628006
		assertEquals(0.005628006, trainer.acceptanceProbability(p), 0.000000001);
		// exp((-0.0000259*200,000,000)/10000) = 0.5957108
		getSettings().setTemperature(10000);
		p = getSettings().getDescriptionParameters().clone();
		assertEquals(0.5957108, trainer.acceptanceProbability(p), 0.000001);
	}

	@Test
	public void testDiffEvalScoreToCurrentlyAcceptedParams() {
		getSettings().setAvgDescriptionEvaluationScore(0.5);
		Parameters p = getSettings().getDescriptionParameters().clone();
		// test first iteration, when accepted Settings are null:
		assertEquals(0.0, trainer.diffEvalScoreToCurrentlyAcceptedParams(p), 0.0);
		// test current Settings equal well performing than accepted ones:
		trainer.setAcceptedParameters(p);
		assertEquals(0.0, trainer.diffEvalScoreToCurrentlyAcceptedParams(p), 0.0);
		// test current Settings better than accepted:
		trainer.setAcceptedParameters(p);
		getSettings().setAvgDescriptionEvaluationScore(1.0);
		p = getSettings().getDescriptionParameters().clone();
		assertEquals(0.5, trainer.diffEvalScoreToCurrentlyAcceptedParams(p), 0.0);
		// test current Settings worse than accepted ones:
		trainer.setAcceptedParameters(p);
		trainer.getAcceptedParameters().setAvgEvaluationScore(0.5);
		getSettings().setAvgDescriptionEvaluationScore(0.25);
		p = getSettings().getDescriptionParameters().clone();
		// 0.25 - 0.5 = -0.25
		assertEquals(-0.25, trainer.diffEvalScoreToCurrentlyAcceptedParams(p),
				0.0);
	}

	@Test
	public void testCoolDown() {
		int t = getSettings().getTemperature();
		trainer.coolDown();
		assertEquals("After cooling down the current temperature should be "
				+ getSettings().getCoolDownBy() + " degrees less than " + t,
				Integer.valueOf(t - getSettings().getCoolDownBy()), getSettings().getTemperature());
	}

	@Test
	public void testAcceptOrRejectParameters() {
		getSettings().setAvgDescriptionEvaluationScore(0.5);
		Parameters p = getSettings().getDescriptionParameters().clone();
		int a = this.trainer.acceptOrRejectParameters(p);
		assertEquals(p,
				this.trainer.getAcceptedParameters());
		assertEquals(
				"The currently evaluated Settings were the first and thus must have been accepted with probability 1.0. Returned int should thus be 3.",
				3, a);
		getSettings().setDescriptionParameters((DescriptionParameters) this.trainer.getAcceptedParameters().neighbour(this.trainer.diffEvalScoreToCurrentlyAcceptedParams(p)));
		getSettings().setAvgDescriptionEvaluationScore(0.75);
		assertTrue(
				"Before calling trainer.acceptOrRejectParameters() accepted Parameters should NOT equal currently evaluated set of Parameters.",
				!getSettings().getDescriptionParameters().equals(
						this.trainer.getAcceptedParameters()));
		p = getSettings().getDescriptionParameters().clone();
		a = this.trainer.acceptOrRejectParameters(p);
		assertEquals(
				"The currently evaluated Settings were better than the currently accepted Settings and thus must have been accepted with probability 1.0. Returned int should thus be 3.",
				3, a);
		assertEquals(getSettings().getDescriptionParameters(),
				this.trainer.getAcceptedParameters());
		// Verify, that some worse performing settings get sometimes accepted or
		// rejected, respectively!
		// P('Accept worse performing Settings') is higher for high Temperatures
		getSettings().setTemperature(10000);
		getSettings().setOptimizationAcceptanceProbabilityScalingFactor(Double.valueOf(1500000));
		p = getSettings().getDescriptionParameters().clone();
		Set<Integer> as = new HashSet<>();
		for (int i = 0; i < 50; i++) {
			this.trainer.getAcceptedParameters().setAvgEvaluationScore(0.75);
			getSettings().setDescriptionParameters((DescriptionParameters) this.trainer.getAcceptedParameters().neighbour(this.trainer.diffEvalScoreToCurrentlyAcceptedParams(p)));
			getSettings().setAvgDescriptionEvaluationScore(0.74538);
			assertTrue(
					"Before calling trainer.acceptOrRejectParameters() accepted Parameters should NOT equal currently evaluated set of Parameters.",
					!getSettings().getDescriptionParameters().equals(this.trainer.getAcceptedParameters()));
			// Verify that we are dealing with the expected
			// acceptance-probability: exp(-0.00462*1500000/10000) = 0.5000736
			p = getSettings().getDescriptionParameters().clone();
			assertEquals(0.5, this.trainer.acceptanceProbability(p), 0.0001);
			as.add(this.trainer.acceptOrRejectParameters(p));
		}
		// Assert, that we did never except a worse parameter set with
		// probability 1.0, which is only applied to better performing
		// parameter-sets!
		String errMsg = "The probability of having only rejections or only acceptances in the"
				+ "above 50 tries can be calculated with the binomial distribution"
				+ "B(0|p=0.5,n=50) < 0.0000000000000009";
		assertTrue(
				"A worse parameter set should have never be accpeted with probability 1.0, which is only applied to better performing parameter-sets",
				!as.contains(3) && !as.contains(2));
		// 50 iterations should have accepted or rejected at least once
		assertTrue(
				"50 iterations should have accepted worse performing settings least once. - "
						+ errMsg, as.contains(1));
		assertTrue(
				"50 iterations should have rejected worse performing settings least once - "
						+ errMsg, as.contains(0));
	}

	/**
	 * Assure that simulated annealing started with the right settings
	 * calculates the average evaluation-score as the pure weighted harmonic
	 * mean of precision and recall, not as difference to best competitors.
	 * 
	 * @throws ParsingException
	 * @throws SAXException
	 * @throws MissingProteinException
	 * @throws MissingAccessionException
	 * @throws IOException
	 */
	@Test
	public void testEvalScoreWithNoCompetitors() throws IOException, MissingAccessionException, MissingProteinException {
		// Default should be FALSE
		assertTrue(!getSettings().getWriteBestBlastHitsToOutput());
		// After setup the competitor annotations
		// "description from best Blast-Hits" should NOT have been saved:
		this.trainer.setup();
		for (Protein p : this.trainer.getProteins().values()) {
			assertTrue(
					"Setting up SimulatedAnnealingTrainer with flag 'write_best_blast_hits_to_output: false' should enforce NOT remembering competitor annotations 'description from best Blast-Hits'.",
					p.getEvaluationScoreCalculator().getBestUnchangedBlastResults()
							.size() == 0);
		}
	}
}
