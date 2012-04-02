package ahrd.controller;

import static ahrd.controller.Settings.getSettings;

import java.io.IOException;
import java.util.HashSet;
import java.util.Random;
import java.util.Set;

import ahrd.exception.MissingInterproResultException;
import ahrd.model.EvaluationScoreCalculator;
import ahrd.model.Protein;
import ahrd.view.TrainerOutputWriter;

public class Trainer extends Evaluator {

	private Parameters acceptedParameters;
	private Parameters bestParameters;
	private Integer bestParametersFoundAtTemperature;
	private TrainerOutputWriter outWriter;
	private Set<Parameters> testedParameters;
	/**
	 * The average of AHRD's maximum evaluation score for each Protein. This is
	 * the maximum of the evaluation scores calculated for all Descriptions of
	 * each Protein. These maximums are then averaged.
	 */
	private Double avgMaxEvaluationScore = 0.0;

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		System.out
				.println("Usage:\njava -Xmx2g -cp ahrd.jar ahrd.controller.Trainer input.yml\n");

		try {
			Trainer trainer = new Trainer(args[0]);
			trainer.setup(false); // false -> Don't log memory and time-usages
			trainer.setupReferences();
			// Blast2GO is another competitor in the field of annotation of
			// predicted Proteins. AHRD might be compared with B2Gs performance:
			trainer.setupBlast2GoAnnots();

			// Try to find optimal parameters heuristically:
			trainer.train();
			// Calculate the average maximum evaluation score AHRD could have
			// possible achieved:
			trainer.calcAvgMaxEvaluationScore();

			// Write final output
			Settings bestSettings = getSettings().clone();
			bestSettings.setParameters(trainer.getBestParameters());
			trainer.outWriter.writeFinalOutput(bestSettings,
					trainer.getAvgMaxEvaluationScore(),
					trainer.getBestParametersFoundAtTemperature());
			System.out
					.println("Logged path through parameter- and score-space into:\n"
							+ getSettings()
									.getPathToSimulatedAnnealingPathLog());
			System.out.println("Written output into:\n"
					+ getSettings().getPathToOutput());
		} catch (Exception e) {
			System.err.println("We are sorry, an unexpected ERROR occurred:");
			e.printStackTrace(System.err);
		}

	}

	/**
	 * Constructor initializes the Settings as given in the argument input.yml
	 * 
	 * @param pathToInputYml
	 * @throws IOException
	 */
	public Trainer(String pathToInputYml) throws IOException {
		super(pathToInputYml);
		this.outWriter = new TrainerOutputWriter();
		// Remember tested Parameter-Sets and their scores?
		if (getSettings().rememberSimulatedAnnealingPath())
			this.testedParameters = new HashSet<Parameters>();
	}

	/**
	 * As of now performs hill-climbing to optimize parameters.
	 * 
	 * @throws IOException
	 * @throws MissingInterproResultException
	 */
	public void train() throws MissingInterproResultException, IOException {
		while (getSettings().getTemperature() > 0) {
			// If we run simulated annealing remembering tested Parameters and
			// their scores,
			// do not calculate current Parameter's performance, if already done
			// in former cycle:
			if (getSettings().rememberSimulatedAnnealingPath()
					&& getTestedParameters().contains(
							getSettings().getParameters())) {
				getSettings().setParameters(
						getAlreadyTestedParameters(getSettings()
								.getParameters()));
			} else {
				// Iterate over all Proteins and assign the best scoring Human
				// Readable Description
				assignHumanReadableDescriptions();
				// Evaluate AHRD's performance for each Protein:
				calculateEvaluationScores();
				// Estimate average performance of current Parameters:
				calcAveragesOfEvalScoreTPRandFPR();
			}
			// Breaking a little bit with the pure simulated annealing
			// algorithm, we remember the best performing Parameters:
			findBestSettings();
			// If started with this option, remember currently evaluated
			// Parameters:
			if (getSettings().rememberSimulatedAnnealingPath())
				getTestedParameters()
						.add(getSettings().getParameters().clone());
			// Initialize the next iteration.
			// Find locally optimal (according to objective function)
			// Parameters:
			int acceptedCurrParameters = acceptOrRejectParameters();
			// Write output of current iteration:
			this.outWriter.writeIterationOutput(getSettings(),
					acceptedCurrParameters);
			// Try a slightly changes set of Parameters:
			initNeighbouringSettings();
			// Cool down temperature:
			coolDown();
		}
	}

	/**
	 * Each iteration the average evaluation-score is compared with the latest
	 * far high-score. If the current Settings Score is better, it will become
	 * the high-score.
	 */
	public void findBestSettings() {
		if (getBestParameters() == null
				|| getSettings().getAvgEvaluationScore() > getBestParameters()
						.getAvgEvaluationScore()) {
			setBestParameters(getSettings().getParameters().clone());
			setBestParametersFoundAtTemperature(getSettings().getTemperature());
		}
	}

	/**
	 * Generates new Settings from the currently accepted ones by
	 * <em>slightly</em> changing them to a <em>neighboring</em> according to
	 * the euclidean distance in the parameter-space Instance.
	 */
	public void initNeighbouringSettings() {
		getSettings().setParameters(getAcceptedParameters().neighbour());
	}

	/**
	 * Calculates Acceptance-Probability according to the <strong>simulated
	 * annealing</strong> algorithm. The distribution of P('Accept worse
	 * performing parameter-sets') := exp(- delta(scores)*scaling-factor /
	 * current-temperature)
	 * 
	 * @return Double - The calculated acceptance-probability
	 */
	public Double acceptanceProbability() {
		// Scaling-Factor referenced for reading convenience. ;-)
		Double sf = getSettings()
				.getOptimizationAcceptanceProbabilityScalingFactor();
		// If current Settings perform better than the so far found best, accept
		// them:
		double p = 1.0;
		// If not, generate Acceptance-Probability based on Score-Difference and
		// current Temperature:
		if (getAcceptedParameters() != null
				&& getSettings().getAvgEvaluationScore() < getAcceptedParameters()
						.getAvgEvaluationScore()) {
			double scoreDiff = getAcceptedParameters().getAvgEvaluationScore()
					- getSettings().getAvgEvaluationScore();
			p = Math.exp(-(scoreDiff * sf) / getSettings().getTemperature());
		}
		return p;
	}

	/**
	 * Diminishes the temperature by one iteration-step.
	 * 
	 * @Note: Temperature is a global Setting.
	 */
	public void coolDown() {
		getSettings().setTemperature(
				getSettings().getTemperature() - getSettings().getCoolDownBy());
	}

	/**
	 * Calculates the average of the difference between AHRD's EvaluationScore
	 * and the best competitor's (objective-function). Also calculates the
	 * average True-Positives- and False-Positives-Rates.
	 */
	public void calcAveragesOfEvalScoreTPRandFPR() {
		// average evaluation-score
		Double avgEvlScr = 0.0;
		// average TPR:
		Double avgTruePosRate = 0.0;
		// average FPR:
		Double avgFalsePosRate = 0.0;
		for (Protein p : getProteins().values()) {
			EvaluationScoreCalculator e = p.getEvaluationScoreCalculator();
			if (e != null) {
				if (e.getEvalScoreMinBestCompScore() != null)
					avgEvlScr += e.getEvalScoreMinBestCompScore();
				if (e.getTruePositivesRate() != null)
					avgTruePosRate += e.getTruePositivesRate();
				if (e.getFalsePositivesRate() != null)
					avgFalsePosRate += e.getFalsePositivesRate();
			}
		}
		// average each number:
		Double numberOfProts = new Double(getProteins().size());
		if (avgEvlScr > 0.0)
			avgEvlScr = avgEvlScr / numberOfProts;
		if (avgTruePosRate > 0.0)
			avgTruePosRate = avgTruePosRate / numberOfProts;
		if (avgFalsePosRate > 0.0)
			avgFalsePosRate = avgFalsePosRate / numberOfProts;
		// done:
		getSettings().setAvgEvaluationScore(avgEvlScr);
		getSettings().setAvgTruePositivesRate(avgTruePosRate);
		getSettings().setAvgFalsePositivesRate(avgFalsePosRate);
	}

	/**
	 * Evaluates the current runs average score (objective function) and based
	 * on this decides according to the simulated annealing algorithm, if the
	 * currently used Settings are accepted or rejected.
	 * 
	 * @Note: Settings are cloned to avoid changing parameters, we want to
	 *        remember unchanged!
	 * 
	 * @return int -
	 *         <ul>
	 *         <li>0 Rejected worse performing parameters</li>
	 *         <li>1 Accepted worse performing parameters</li>
	 *         <li>2 Accepted equally well performing parameters</li>
	 *         <li>3 Accepted better performing parameters</li>
	 *         </ul>
	 */
	public int acceptOrRejectParameters() {
		int accepted = 0; // Rejected worse performing parameters
		double acceptCurrSettingsProb = acceptanceProbability();
		if (acceptCurrSettingsProb == 1.0) {
			if (getAcceptedParameters() == null
					|| getAcceptedParameters().getAvgEvaluationScore() < getSettings()
							.getAvgEvaluationScore()) {
				accepted = 3; // Accepted better performing parameters
			} else {
				accepted = 2; // Accepted equally well performing parameters
			}
			setAcceptedParameters(getSettings().getParameters().clone());
		} else {
			// Take random decision
			Random r = Utils.random;
			if (r.nextDouble() <= acceptCurrSettingsProb) {
				setAcceptedParameters(getSettings().getParameters().clone());
				accepted = 1; // Accepted worse performing parameters
			}
			// else discard the current Settings and continue with the so far
			// optimal ones.
		}
		return accepted;
	}

	/**
	 * Do not calculate the current Parameters' performance again, use
	 * remembered scores instead.
	 * 
	 * @param current
	 * @return Parameters
	 */
	public Parameters getAlreadyTestedParameters(Parameters current) {
		Parameters alreadyTested = null;
		for (Parameters iterParams : getTestedParameters().toArray(
				new Parameters[] {})) {
			if (current.equals(iterParams))
				alreadyTested = iterParams;
		}
		return alreadyTested;
	}

	/**
	 * This calculates the average maximum evaluation score AHRD could possibly
	 * achieve.
	 */
	public void calcAvgMaxEvaluationScore() {
		for (Protein prot : getProteins().values()) {
			prot.getEvaluationScoreCalculator()
					.findHighestPossibleEvaluationScore();
			setAvgMaxEvaluationScore(getAvgMaxEvaluationScore()
					+ prot.getEvaluationScoreCalculator()
							.getHighestPossibleEvaluationScore());
		}
		setAvgMaxEvaluationScore(getAvgMaxEvaluationScore()
				/ getProteins().size());
	}

	public Parameters getAcceptedParameters() {
		return acceptedParameters;
	}

	public void setAcceptedParameters(Parameters acceptedSettings) {
		this.acceptedParameters = acceptedSettings;
	}

	public Parameters getBestParameters() {
		return bestParameters;
	}

	public void setBestParameters(Parameters bestParameters) {
		this.bestParameters = bestParameters;
	}

	public Set<Parameters> getTestedParameters() {
		return testedParameters;
	}

	public Double getAvgMaxEvaluationScore() {
		return avgMaxEvaluationScore;
	}

	public void setAvgMaxEvaluationScore(Double avgMaxEvaluationScore) {
		this.avgMaxEvaluationScore = avgMaxEvaluationScore;
	}

	public Integer getBestParametersFoundAtTemperature() {
		return bestParametersFoundAtTemperature;
	}

	public void setBestParametersFoundAtTemperature(
			Integer bestParametersFoundAtTemperature) {
		this.bestParametersFoundAtTemperature = bestParametersFoundAtTemperature;
	}

}
