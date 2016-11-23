package ahrd.test;

import static ahrd.controller.Settings.getSettings;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.sql.SQLException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

import nu.xom.ParsingException;

import org.junit.Before;
import org.junit.Test;
import org.xml.sax.SAXException;

import ahrd.controller.Parameters;
import ahrd.controller.Trainer;
import ahrd.exception.MissingAccessionException;
import ahrd.exception.MissingInterproResultException;
import ahrd.exception.MissingProteinException;
import ahrd.model.Protein;

public class TrainerTest {

	private Trainer trainer;

	@Before
	public void setUp() throws IOException, MissingAccessionException,
			MissingProteinException, SAXException, ParsingException {
		trainer = new Trainer("./test/resources/trainer_input.yml");
		trainer.setup(false); // false -> Don't log memory and time-usages
		trainer.setupReferenceDescriptions();
		// Sets up description and GOterm annotations of competitors in the EvaluationScoreCalculators of each Protein
		trainer.setupCompetitors();
	}

	@Test
	public void testAvgEvaluationScore() {
		assertTrue("Trainer should initialize Settings to Training-Mode.",
				getSettings().isInTrainingMode());
		Protein p1 = new Protein("protein_one", 200);
		p1.getEvaluationScoreCalculator().setEvalutionScore(1.0);
		p1.getEvaluationScoreCalculator().setEvalScoreMinBestCompScore(1.0);
		p1.getEvaluationScoreCalculator().setTruePositivesRate(0.6);
		p1.getEvaluationScoreCalculator().setFalsePositivesRate(0.7);
		Protein p2 = new Protein("protein_two", 210);
		p2.getEvaluationScoreCalculator().setEvalutionScore(0.8);
		p2.getEvaluationScoreCalculator().setEvalScoreMinBestCompScore(0.8);
		p2.getEvaluationScoreCalculator().setTruePositivesRate(0.7);
		p2.getEvaluationScoreCalculator().setFalsePositivesRate(0.5);
		Protein p3 = new Protein("protein_three", 220);
		p3.getEvaluationScoreCalculator().setEvalutionScore(0.3);
		p3.getEvaluationScoreCalculator().setEvalScoreMinBestCompScore(0.3);
		p3.getEvaluationScoreCalculator().setTruePositivesRate(0.5);
		p3.getEvaluationScoreCalculator().setFalsePositivesRate(0.3);
		this.trainer.setProteins(new HashMap<String, Protein>());
		this.trainer.getProteins().put(p1.getAccession(), p1);
		this.trainer.getProteins().put(p2.getAccession(), p2);
		this.trainer.getProteins().put(p3.getAccession(), p3);
		// test
		trainer.calcAveragesOfEvalScoreTPRandFPR();
		// (1.0 + 0.8 + 0.3) / 3 = 0.7
		assertEquals(0.7, getSettings().getAvgEvaluationScore(), 0.000000000001);
		// (0.6 + 0.7 + 0.5) / 3 = 0.6
		assertEquals(0.6, getSettings().getAvgTruePositivesRate(),
				0.000000000001);
		// (0.7 + 0.5 + 0.3) / 3 = 0.5
		assertEquals(0.5, getSettings().getAvgFalsePositivesRate(),
				0.000000000001);
		// test robustness for zero scores:
		for (Protein p : trainer.getProteins().values()) {
			p.getEvaluationScoreCalculator().setEvalScoreMinBestCompScore(0.0);
		}
		trainer.calcAveragesOfEvalScoreTPRandFPR();
		assertEquals(getSettings().getAvgEvaluationScore(), 0.7, 0.00000000001);
	}

	@Test
	public void testRememberSimulatedAnnealingPath()
			throws MissingInterproResultException, IOException, SQLException {
		// Just do two cycles:
		getSettings().setTemperature(2);
		Parameters p = getSettings().getParameters().clone();
		this.trainer.train();
		// test:
		assertEquals(
				"Training should have visited two distinct parameter-sets.", 2,
				this.trainer.getTestedParameters().size());
		Parameters clone = this.trainer.getAlreadyTestedParameters(p);
		assertNotNull(clone);
		// Scores should be equal:
		assertNotNull("Avg EvaluationScore should be remembered.",
				clone.getAvgEvaluationScore());
		assertNotNull("Avg TPR should be remembered.",
				clone.getAvgTruePositivesRate());
		assertNotNull("Avg FPR should be remembered.",
				clone.getAvgFalsePositivesRate());
	}

	@Test
	public void testAcceptanceProbability() {
		getSettings().setAvgEvaluationScore(0.5);
		getSettings().setOptimizationAcceptanceProbabilityScalingFactor(
				new Double(200000000));
		getSettings().setTemperature(1000);
		// test first iteration, when accepted Settings are null:
		assertEquals(1.0, trainer.acceptanceProbability(), 0.0);
		// test current Settings better than accepted:
		trainer.setAcceptedParameters(getSettings().getParameters().clone());
		getSettings().setAvgEvaluationScore(1.0);
		assertEquals(1.0, trainer.acceptanceProbability(), 0.0);
		// test current Settings worse than accepted ones:
		trainer.setAcceptedParameters(getSettings().getParameters().clone());
		getSettings().setAvgEvaluationScore(0.9999741);
		// 0.9999741 - 1.0 = -2.59 * 10^-5
		assertEquals(-0.0000259,
				trainer.diffEvalScoreToCurrentlyAcceptedParams(), 0.00000001);
		// exp((-0.0000259*200,000,000)/1000) = 0.005628006
		assertEquals(0.005628006, trainer.acceptanceProbability(), 0.000000001);
		// exp((-0.0000259*200,000,000)/10000) = 0.5957108
		getSettings().setTemperature(10000);
		assertEquals(0.5957108, trainer.acceptanceProbability(), 0.000001);
	}

	@Test
	public void testDiffEvalScoreToCurrentlyAcceptedParams() {
		getSettings().setAvgEvaluationScore(0.5);
		// test first iteration, when accepted Settings are null:
		assertEquals(0.0, trainer.diffEvalScoreToCurrentlyAcceptedParams(), 0.0);
		// test current Settings equal well performing than accepted ones:
		trainer.setAcceptedParameters(getSettings().getParameters().clone());
		assertEquals(0.0, trainer.diffEvalScoreToCurrentlyAcceptedParams(), 0.0);
		// test current Settings better than accepted:
		trainer.setAcceptedParameters(getSettings().getParameters().clone());
		getSettings().setAvgEvaluationScore(1.0);
		assertEquals(0.5, trainer.diffEvalScoreToCurrentlyAcceptedParams(), 0.0);
		// test current Settings worse than accepted ones:
		trainer.setAcceptedParameters(getSettings().getParameters().clone());
		trainer.getAcceptedParameters().setAvgEvaluationScore(0.5);
		getSettings().setAvgEvaluationScore(0.25);
		// 0.25 - 0.5 = -0.25
		assertEquals(-0.25, trainer.diffEvalScoreToCurrentlyAcceptedParams(),
				0.0);
	}

	@Test
	public void testCoolDown() {
		Integer t = new Integer(getSettings().getTemperature());
		trainer.coolDown();
		assertEquals("After cooling down the current temperature should be "
				+ getSettings().getCoolDownBy() + " degrees less than " + t,
				new Integer(t - getSettings().getCoolDownBy()), getSettings()
						.getTemperature());
	}

	@Test
	public void testAcceptOrRejectParameters() {
		getSettings().setAvgEvaluationScore(0.5);
		int a = this.trainer.acceptOrRejectParameters();
		assertEquals(getSettings().getParameters(),
				this.trainer.getAcceptedParameters());
		assertEquals(
				"The currently evaluated Settings were the first and thus must have been accepted with probability 1.0. Returned int should thus be 3.",
				3, a);
		this.trainer.initNeighbouringSettings();
		getSettings().setAvgEvaluationScore(0.75);
		assertTrue(
				"Before calling trainer.acceptOrRejectParameters() accepted Parameters should NOT equal currently evaluated set of Parameters.",
				!getSettings().getParameters().equals(
						this.trainer.getAcceptedParameters()));
		a = this.trainer.acceptOrRejectParameters();
		assertEquals(
				"The currently evaluated Settings were better than the currently accepted Settings and thus must have been accepted with probability 1.0. Returned int should thus be 3.",
				3, a);
		assertEquals(getSettings().getParameters(),
				this.trainer.getAcceptedParameters());
		// Verify, that some worse performing settings get sometimes accepted or
		// rejected, respectively!
		// P('Accept worse performing Settings') is higher for high Temperatures
		getSettings().setTemperature(10000);
		getSettings().setOptimizationAcceptanceProbabilityScalingFactor(
				new Double(1500000));
		Set<Integer> as = new HashSet<Integer>();
		for (int i = 0; i < 50; i++) {
			this.trainer.getAcceptedParameters().setAvgEvaluationScore(0.75);
			this.trainer.initNeighbouringSettings();
			getSettings().setAvgEvaluationScore(0.74538);
			assertTrue(
					"Before calling trainer.acceptOrRejectParameters() accepted Parameters should NOT equal currently evaluated set of Parameters.",
					!getSettings().getParameters().equals(
							this.trainer.getAcceptedParameters()));
			// Verify that we are dealing with the expected
			// acceptance-probability: exp(-0.00462*1500000/10000) = 0.5000736
			assertEquals(0.5, this.trainer.acceptanceProbability(), 0.0001);
			as.add(this.trainer.acceptOrRejectParameters());
		}
		// Assert, that we did never except a worse parameter set with
		// probability 1.0, which is only applied to better performing
		// parameter-sets!
		String errMsg = "The probability of having only rejections or only acceptances in the"
				+ "above 50 tries can be calculated with the binomial distribution"
				+ "B(0|p=0.5,n=50) < 0.0000000000000009";
		assertTrue(
				"A worse parameter set should have never be accpeted with probability 1.0, which is only applied to better performing parameter-sets",
				!as.contains(3.0) && !as.contains(2.0));
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
	public void testEvalScoreWithNoCompetitors() throws IOException,
			MissingAccessionException, MissingProteinException, SAXException,
			ParsingException {
		// Default should be FALSE
		assertTrue(!getSettings().getWriteBestBlastHitsToOutput());
		// After setup the competitor annotations
		// "description from best Blast-Hits" should NOT have been saved:
		this.trainer.setup(false);
		for (Protein p : this.trainer.getProteins().values()) {
			assertTrue(
					"Setting up Trainer with flag 'write_best_blast_hits_to_output: false' should enforce NOT remembering competitor annotations 'description from best Blast-Hits'.",
					p.getEvaluationScoreCalculator().getBestUnchangedBlastResults()
							.size() == 0);
		}
	}
}
