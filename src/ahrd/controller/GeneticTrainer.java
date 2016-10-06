package ahrd.controller;

import static ahrd.controller.Settings.getSettings;

import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.NavigableSet;
import java.util.Random;
import java.util.Set;
import java.util.TreeSet;

import ahrd.exception.MissingInterproResultException;
import ahrd.model.EvaluationScoreCalculator;
import ahrd.model.Protein;
import ahrd.view.TrainerOutputWriter;

public class GeneticTrainer extends Evaluator {

	private static final int NUMBER_OF_GENERATIONS = 5;
	private static final int POPULATION_SIZE = 10;
	private static final Double GENERATIONAL_SURVIVAL_RATE = 0.2;
	private static final Double GENERATIONAL_OFFSPRING_RATE = 0.2;
	private static final Double GENERATIONAL_MUTANT_RATE = 0.2;

	private static int numberOfSurvivors = (int) Math.round(POPULATION_SIZE * GENERATIONAL_SURVIVAL_RATE);
	private static int numberOfOffspring = (int) Math.round(POPULATION_SIZE * GENERATIONAL_OFFSPRING_RATE);
	private static int numberOfMutants = (int) Math.round(POPULATION_SIZE * GENERATIONAL_MUTANT_RATE);
	private Parameters bestParameters;
	private Integer generationBestParametersWereFoundIn;
	private TrainerOutputWriter outWriter;
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
		System.out.println("Usage:\njava -Xmx2g -cp ahrd.jar ahrd.controller.GeneticTrainer input.yml\n");

		try {
			GeneticTrainer trainer = new GeneticTrainer(args[0]);
			trainer.setup(false); // false -> Don't log memory and time-usages
			// After the setup the unique short accessions are no longer needed:
			trainer.setUniqueBlastResultShortAccessions(null);
			trainer.setupReferenceDescriptions();

			// Try to find optimal parameters heuristically:
			trainer.train();
			// Calculate the average maximum evaluation score AHRD could have
			// possible achieved:
			trainer.calcAvgMaxEvaluationScore();

			// Write final output
			Settings bestSettings = getSettings().clone();
			bestSettings.setParameters(trainer.getBestParameters());
			trainer.outWriter.writeFinalOutput(bestSettings, trainer.getAvgMaxEvaluationScore(),
					trainer.getGenerationBestParametersWereFoundIn());
			System.out.println("Logged path through parameter- and score-space into:\n"
					+ getSettings().getPathToSimulatedAnnealingPathLog());
			System.out.println("Written output into:\n" + getSettings().getPathToOutput());
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
	public GeneticTrainer(String pathToInputYml) throws IOException {
		super(pathToInputYml);
		this.outWriter = new TrainerOutputWriter();
		// Remember tested Parameter-Sets and their scores?
	}

	/**
	 * As of now performs hill-climbing to optimize parameters.
	 * 
	 * @throws IOException
	 * @throws MissingInterproResultException
	 * @throws SQLException
	 */
	public void train() throws MissingInterproResultException, IOException, SQLException {
		// Set up first generation
		Set<Parameters> population = new HashSet<Parameters>();
		for (int i = 1; i <= POPULATION_SIZE; i++) {
			List<String> sortedDistinctBlastDatabaseNames = new ArrayList<String>();
			sortedDistinctBlastDatabaseNames.addAll(getSettings().getBlastDatabases());
			Collections.sort(sortedDistinctBlastDatabaseNames);
			population.add(Parameters.randomParameters(sortedDistinctBlastDatabaseNames));
		}
		int generation = 1;
		while (generation <= NUMBER_OF_GENERATIONS) {
			// Determine the fitness of each individual (parameter set) in the
			// population
			for (Iterator<Parameters> paraIter = population.iterator(); paraIter.hasNext();) {
				Parameters individual = paraIter.next();
				getSettings().setParameters(individual);
				// Iterate over all Proteins and assign the best scoring Human
				// Readable Description
				if (individual.getAvgEvaluationScore() == null) {
					assignHumanReadableDescriptions();
					// Evaluate AHRD's performance for each Protein:
					calculateEvaluationScores();
					// Estimate average performance of current Parameters:
					calcAveragesOfEvalScoreTPRandFPR();
				}
			}
			// Survival of the fittest
			NavigableSet<Parameters> ranking = new TreeSet<Parameters>();
			ranking.addAll(population);
			population.clear();
			for (int i = 1; i <= numberOfSurvivors; i++) {
				population.add(ranking.pollLast());
			}
			/*
			 * OR: pollFirst until ranking.size == numberOfSurvivors -> keep
			 * ranking of survivors to simulate sexual selection (bias toward
			 * higher scores when randomly selecting mates) -> keep ranking of
			 * survivors for bias in selection of perm to mutate
			 */
			population.addAll(ranking);
			// Recombination of survivors
			ranking.clear();
			ranking.addAll(population);
			for (int i = 1; i <= numberOfOffspring; i++) {
				Random rand = Utils.random;
				// population.
				// rand.nextInt(numberOfSurvivors)+1;
				// population.add(e);
			}

			// Increment generation counter
			generation += 1;
		}
	}

	/**
	 * Calculates the average of AHRD's EvaluationScore (objective-function).
	 * Also calculates the average True-Positives- and False-Positives-Rates.
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
				if (e.getEvalutionScore() != null)
					avgEvlScr += e.getEvalutionScore();
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
	 * This calculates the average maximum evaluation score AHRD could possibly
	 * achieve.
	 */
	public void calcAvgMaxEvaluationScore() {
		for (Protein prot : getProteins().values()) {
			prot.getEvaluationScoreCalculator().findHighestPossibleEvaluationScore();
			setAvgMaxEvaluationScore(getAvgMaxEvaluationScore()
					+ prot.getEvaluationScoreCalculator().getHighestPossibleEvaluationScore());
		}
		setAvgMaxEvaluationScore(getAvgMaxEvaluationScore() / getProteins().size());
	}

	public Double getAvgMaxEvaluationScore() {
		return avgMaxEvaluationScore;
	}

	public void setAvgMaxEvaluationScore(Double avgMaxEvaluationScore) {
		this.avgMaxEvaluationScore = avgMaxEvaluationScore;
	}

	public Parameters getBestParameters() {
		return bestParameters;
	}

	public void setBestParameters(Parameters bestParameters) {
		this.bestParameters = bestParameters;
	}

	public Integer getGenerationBestParametersWereFoundIn() {
		return generationBestParametersWereFoundIn;
	}

	public void setGenerationBestParametersWereFoundIn(Integer generationBestParametersWereFoundIn) {
		this.generationBestParametersWereFoundIn = generationBestParametersWereFoundIn;
	}

}
