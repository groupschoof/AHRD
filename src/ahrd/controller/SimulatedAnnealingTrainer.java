package ahrd.controller;

import static ahrd.controller.Settings.getSettings;

import java.io.IOException;
import java.sql.SQLException;
import java.util.HashSet;
import java.util.Random;
import java.util.Set;

import ahrd.exception.MissingAccessionException;
import ahrd.exception.MissingInterproResultException;
import ahrd.view.TrainerOutputWriter;

public class SimulatedAnnealingTrainer extends Trainer {

	private Parameters acceptedParameters;
	private Integer bestParametersFoundAtTemperature;
	private Set<Parameters> testedParameters;

	/**
	 * Constructor initializes the Settings as given in the argument input.yml
	 * 
	 * @param pathToInputYml
	 * @throws IOException
	 */
	public SimulatedAnnealingTrainer(String pathToInputYml) throws IOException {
		super(pathToInputYml);
		this.outWriter = new TrainerOutputWriter();
		// Remember tested Parameter-Sets and their scores?
		if (getSettings().rememberSimulatedAnnealingPath())
			this.testedParameters = new HashSet<Parameters>();
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		System.out
				.println("Usage:\njava -Xmx2g -cp ahrd.jar ahrd.controller.SimulatedAnnealingTrainer input.yml\n");

		try {
			SimulatedAnnealingTrainer trainer = new SimulatedAnnealingTrainer(args[0]);
			trainer.setup(false); // false -> Don't log memory and time-usages
			if (getSettings().hasGeneOntologyAnnotations() && getSettings().hasGroundTruthGoAnnotations()) {
				getSettings().setFindHighestPossibleGoScore(true);
			}
			// After the setup the unique short accessions are no longer needed:
			trainer.setUniqueBlastResultShortAccessions(null);
			trainer.setupGroundTruthDescriptions();
			trainer.setupGoAnnotationEvaluation();
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
	 * As of now performs hill-climbing to optimize parameters.
	 * 
	 * @throws IOException
	 * @throws MissingInterproResultException
	 * @throws SQLException
	 * @throws MissingAccessionException 
	 */
	public void train() throws MissingInterproResultException, IOException,
			SQLException, MissingAccessionException {
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
				reinitializeBlastResults();
				// Iterate over all Proteins and assign the best scoring Human
				// Readable Description
				assignHumanReadableDescriptions();
				if (getSettings().hasGeneOntologyAnnotations() && getSettings().hasGroundTruthGoAnnotations()) {
					assignGeneOntologyTerms();
					goAnnotsStringToObject();
				}
				// Evaluate AHRD's performance for each Protein:
				calculateEvaluationScores();
				// Estimate average performance of current Parameters:
				calcAveragesOfEvalScorePrecisionAndRecall();
			}
			// Breaking a little bit with the pure simulated annealing
			// algorithm, we remember the best performing Parameters:
			findBestSettings();
			// If started with this option, remember currently evaluated
			// Parameters:
			if (getSettings().rememberSimulatedAnnealingPath())
				getTestedParameters()
						.add(getSettings().getParameters().clone());
			// Remember difference in avg. evaluation-scores, *before* accepting
			// or rejecting current Parameters:
			Double diffScores = diffEvalScoreToCurrentlyAcceptedParams();
			// Initialize the next iteration.
			// Find locally optimal (according to objective function)
			// Parameters:
			int acceptedCurrParameters = acceptOrRejectParameters();
			// Write output of current iteration:
			this.outWriter.writeIterationOutput(getSettings(), diffScores,
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
		getSettings().setParameters(
				getAcceptedParameters().neighbour(
						diffEvalScoreToCurrentlyAcceptedParams()));
	}

	public Double diffEvalScoreToCurrentlyAcceptedParams() {
		return (getAcceptedParameters() != null) ? getSettings()
				.getAvgEvaluationScore()
				- getAcceptedParameters().getAvgEvaluationScore() : 0.0;
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
				&& diffEvalScoreToCurrentlyAcceptedParams() < 0.0) {
			// In this case the difference in avg. evaluation scores of current
			// to accepted parameters is always NEGATIVE.
			// Hence the following formula can be written as:
			// p := exp((delta.scores*sf)/T.curr), where delta.score is a
			// negative real value.
			p = Math.exp(diffEvalScoreToCurrentlyAcceptedParams() * sf
					/ getSettings().getTemperature());
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

	public Parameters getAcceptedParameters() {
		return acceptedParameters;
	}

	public void setAcceptedParameters(Parameters acceptedSettings) {
		this.acceptedParameters = acceptedSettings;
	}

	public Set<Parameters> getTestedParameters() {
		return testedParameters;
	}

	public Integer getBestParametersFoundAtTemperature() {
		return bestParametersFoundAtTemperature;
	}

	public void setBestParametersFoundAtTemperature(
			Integer bestParametersFoundAtTemperature) {
		this.bestParametersFoundAtTemperature = bestParametersFoundAtTemperature;
	}

}
