package ahrd.controller;

import static ahrd.controller.Settings.getSettings;

import java.io.IOException;
import java.util.HashSet;
import java.util.Random;
import java.util.Set;

import org.semanticweb.owlapi.model.OWLOntologyCreationException;

import ahrd.exception.MissingAccessionException;
import ahrd.view.SimulatedAnnealingTrainerOutputWriter;

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
		// Remember tested Parameter-Sets and their scores?
		if (getSettings().rememberSimulatedAnnealingPath())	this.testedParameters = new HashSet<>();
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		System.out
				.println("Usage:\njava -cp ahrd.jar ahrd.controller.SimulatedAnnealingTrainer input.yml\n");
		try {
			// Try to heuristically find optimal parameters for the annotation with descriptions
			SimulatedAnnealingTrainer trainer = new SimulatedAnnealingTrainer(args[0]);
			if (getSettings().doEvaluateDescriptions()) {
				if (getSettings().doAnnotateGoTerms()) {
					getSettings().setAnnotateGoTerms(false); // prevent parsing of GO references (can save a LOT of time)
					trainer.setup(true); // true -> Log memory and time-usages
					getSettings().setAnnotateGoTerms(true);
				} else {
					trainer.setup(true); // true -> Log memory and time-usages
				}
				trainer.setUniqueBlastResultShortAccessions(null); // After the setup the unique short accessions are no longer needed
				trainer.setupGroundTruthDescriptions();
				trainer.outputWriter = new SimulatedAnnealingTrainerOutputWriter(getSettings().getPathToDescriptionTrainingPathLog());
				Parameters seed = getSettings().getDescriptionParameters().clone();
				trainer.outputWriter.writeHeader(seed);
				trainer.train(seed);
				trainer.calcAvgMaxDescriptionScore();
				getSettings().setPathToOutput(getSettings().getPathToDescriptionOutput());
				trainer.outputWriter.writeFinalOutput(
						trainer.getBestParametersFoundAtTemperature(),
						trainer.getAvgMaxDescriptionScore(),
						trainer.getBestParameters());
			}
			// Try to heuristically find optimal parameters for the annotation with GO terms
			trainer = new SimulatedAnnealingTrainer(args[0]);
			if (getSettings().doEvaluateGoTerms()) { 
				trainer.setup(true); // true -> Log memory and time-usages
				trainer.setUniqueBlastResultShortAccessions(null); // After the setup the unique short accessions are no longer needed
				getSettings().setFindHighestPossibleGoScore(true);
				trainer.setupGoAnnotationEvaluation();
				trainer.outputWriter = new SimulatedAnnealingTrainerOutputWriter(getSettings().getPathToGoTrainingPathLog());
				Parameters seed = getSettings().getGoParameters().clone();
				trainer.outputWriter.writeHeader(seed);
				trainer.train(seed);
				trainer.calcAvgMaxGoScore();
				getSettings().setPathToOutput(getSettings().getPathToGoOutput());
				trainer.outputWriter.writeFinalOutput(
						trainer.getBestParametersFoundAtTemperature(),
						trainer.getAvgMaxGoScore(),
						trainer.getBestParameters());
			}
			System.out.println("Logged path through parameter- and score-space into:\n" + getSettings().getPathToDescriptionTrainingPathLog());
			System.out.println("Written output into:\n"	+ getSettings().getPathToOutput());
		} catch (Exception e) {
			System.err.println("We are sorry, an unexpected ERROR occurred:");
			e.printStackTrace(System.err);
		}

	}

	/**
	 * As of now performs hill-climbing to optimize parameters.
	 * @throws IOException 
	 * @throws MissingAccessionException 
	 * @throws OWLOntologyCreationException 
	 * 
	 */
	public void train(Parameters currentParameters) throws IOException, OWLOntologyCreationException, MissingAccessionException {
		while (getSettings().getTemperature() > 0) {
			// If we run simulated annealing remembering tested Parameters and their scores,
			// do not calculate current Parameter's performance, if already done in former cycle.
			// (The new Parameters object has only the parameters.
			// The old one on the other hand also has the evaluation scores. 
			if (getSettings().rememberSimulatedAnnealingPath() && getTestedParameters().contains(currentParameters)) {
				currentParameters = getAlreadyTestedParameters(currentParameters);
			} else {
				if (currentParameters instanceof DescriptionParameters) { // We are optimizing description annotation
					getSettings().setDescriptionParameters((DescriptionParameters) currentParameters);
					reinitializeBlastResults();
					// Iterate over all Proteins and assign the best scoring Human Readable Description
					assignHumanReadableDescriptions();
					// Evaluate AHRD's performance for each Protein:
					calculateEvaluationScores();
					// Estimate average performance of current Parameters:
					calcAveragesOfDescriptionScores();
				} else { // We are optimizing GO annotation
					getSettings().setGoParameters((GoParameters) currentParameters);
					reinitializeBlastResults();
					// Iterate over all Proteins and assign the best scoring GO terms
					assignGeneOntologyTerms();
					goAnnotsStringToObject();
					// Evaluate AHRD's performance for each Protein:
					calculateEvaluationScores();
					// Estimate average performance of current Parameters:
					calcAveragesOfGoScores();
				}
			}
			
			double diff;
			// Breaking a little bit with the pure simulated annealing algorithm, we remember the best performing Parameters
			if (getBestParameters() == null || currentParameters.getAvgEvaluationScore() > getBestParameters().getAvgEvaluationScore()) {
				setBestParameters(currentParameters.clone());
				setBestParametersFoundAtTemperature(getSettings().getTemperature());
			}
			// If started with this option, remember currently evaluated parameters:
			if (getSettings().rememberSimulatedAnnealingPath()) {
				getTestedParameters().add(currentParameters.clone());
			}
			// Remember difference in avg. evaluation-scores, *before* accepting or rejecting current Parameters:
			diff = diffEvalScoreToCurrentlyAcceptedParams(currentParameters);
			// Initialize the next iteration.
			// Find locally optimal (according to objective function) parameters:
			Integer acceptedCurrParameters = acceptOrRejectParameters(currentParameters);
			// Write output of current iteration
			this.outputWriter.writeIterationOutput(
					getSettings().getTemperature(), currentParameters, diff, acceptedCurrParameters.toString());
			// Try a slightly changed set of Parameters:
			currentParameters = getAcceptedParameters().neighbour(diffEvalScoreToCurrentlyAcceptedParams(currentParameters));
			// Cool down temperature:
			coolDown();
		}
	}

	public Double diffEvalScoreToCurrentlyAcceptedParams(Parameters currentParameters) {
		Double diffScores = 0.0;
		if (getAcceptedParameters() != null) {
			diffScores = currentParameters.getAvgEvaluationScore() - getAcceptedParameters().getAvgEvaluationScore();
		}
		return diffScores;
	}
	
	/**
	 * Calculates Acceptance-Probability according to the <strong>simulated
	 * annealing</strong> algorithm. The distribution of P('Accept worse
	 * performing parameter-sets') := exp(- delta(scores)*scaling-factor /
	 * current-temperature)
	 * 
	 * @return Double - The calculated acceptance-probability
	 */
	public Double acceptanceProbability(Parameters curretParameters) {
		// Scaling-Factor referenced for reading convenience. ;-)
		Double sf = getSettings().getOptimizationAcceptanceProbabilityScalingFactor();
		// If current Settings perform better than the so far found best, accept
		// them:
		double p = 1.0;
		// If not, generate Acceptance-Probability based on Score-Difference and
		// current Temperature:
		if (getAcceptedParameters() != null	&& diffEvalScoreToCurrentlyAcceptedParams(curretParameters) < 0.0) {
			// In this case the difference in avg. evaluation scores of current
			// to accepted parameters is always NEGATIVE.
			// Hence the following formula can be written as:
			// p := exp((delta.scores*sf)/T.curr), where delta.score is a
			// negative real value.
			p = Math.exp(diffEvalScoreToCurrentlyAcceptedParams(curretParameters) * sf / getSettings().getTemperature());
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
	 * @throws CloneNotSupportedException 
	 */
	public int acceptOrRejectParameters(Parameters curretParameters) {
		int accepted = 0; // Rejected worse performing parameters
		double acceptCurrSettingsProb = acceptanceProbability(curretParameters);
		if (acceptCurrSettingsProb == 1.0) {
			if (getAcceptedParameters() == null
					|| getAcceptedParameters().getAvgEvaluationScore() < curretParameters.getAvgEvaluationScore()) {
				accepted = 3; // Accepted better performing parameters
			} else {
				accepted = 2; // Accepted equally well performing parameters
			}
			setAcceptedParameters(curretParameters.clone());
		} else {
			// Take random decision
			Random r = Utils.random;
			if (r.nextDouble() <= acceptCurrSettingsProb) {
				setAcceptedParameters(curretParameters.clone());
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

	public Integer getBestParametersFoundAtTemperature() {
		return bestParametersFoundAtTemperature;
	}

	public void setBestParametersFoundAtTemperature(int temperature) {
		this.bestParametersFoundAtTemperature = temperature;
	}
	
	public Set<Parameters> getTestedParameters() {
		return testedParameters;
	}

}
