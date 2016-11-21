package ahrd.controller;

import static ahrd.controller.Settings.getSettings;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.NavigableSet;
import java.util.Set;
import java.util.TreeSet;

import ahrd.exception.MissingAccessionException;
import ahrd.exception.MissingInterproResultException;
import ahrd.exception.MissingProteinException;
import ahrd.model.EvaluationScoreCalculator;
import ahrd.model.Protein;
import ahrd.view.GeneticTrainerOutputWriter;
import nu.xom.ParsingException;

public class GeneticTrainer extends Evaluator {

	private static final Double GENERATIONAL_SURVIVAL_RATE = 0.2;
	private static final Double GENERATIONAL_OFFSPRING_RATE = 0.2;
	private static final Double GENERATIONAL_MUTANT_RATE = 0.2;

	private static int numberOfSurvivors;
	private static int numberOfOffspring;
	private static int numberOfMutants;
	private Parameters bestParameters;
	private Integer generationBestParametersWereFoundIn;
	private GeneticTrainerOutputWriter outWriter;
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
			trainer.setupGoAnnotationEvaluation();
			// Try to find optimal parameters heuristically:
			trainer.train();
			// Calculate the average maximum evaluation score AHRD could have
			// possibly achieved:
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
		this.outWriter = new GeneticTrainerOutputWriter();
		numberOfSurvivors = (int) Math.round(getSettings().getPopulationSize() * GENERATIONAL_SURVIVAL_RATE);
		numberOfOffspring = (int) Math.round(getSettings().getPopulationSize() * GENERATIONAL_OFFSPRING_RATE);
		numberOfMutants = (int) Math.round(getSettings().getPopulationSize() * GENERATIONAL_MUTANT_RATE);
	}

	/**
	 * Random first generation. Survival of the fittest, recombination of fit
	 * survivors, mutants of fit survivors and random parameter sets for the
	 * rest in each succeeding generation.
	 * 
	 * @throws IOException
	 * @throws MissingInterproResultException
	 * @throws SQLException
	 * @throws MissingAccessionException 
	 * @throws ParsingException 
	 * @throws MissingProteinException 
	 */
	public void train() throws MissingInterproResultException, IOException, SQLException, MissingAccessionException, MissingProteinException, ParsingException {
		Set<Parameters> population = new HashSet<Parameters>();
		// Set up first generation
		List<String> sortedDistinctBlastDatabaseNames = new ArrayList<String>();
		sortedDistinctBlastDatabaseNames.addAll(getSettings().getBlastDatabases());
		Collections.sort(sortedDistinctBlastDatabaseNames);
		// Add parameters from YML-input to first generation (enables seeding with high mean evaluation score parameter set)
		Parameters seed = getSettings().getParameters().clone();
		seed.setOrigin("seed");
		population.add(seed);
		for (int i = 2; i <= getSettings().getPopulationSize(); i++) {
			population.add(Parameters.randomParameters(sortedDistinctBlastDatabaseNames));
		}
		int generation = 1;
		double diffAvgEvalScoreToLastGeneration = 0;
		// simulate generational succession
		while (generation <= getSettings().getNumberOfGenerations()) {
			// Show progress
			System.out.print("\rEvaluating generation " + generation + " of " + getSettings().getNumberOfGenerations());
			// Determine the fitness of each individual (parameter set) in the
			// population
			for (Parameters individual : population) {
//				if (individual.getAvgEvaluationScore() == null) {
					getSettings().setParameters(individual);
					// Iterate over all Proteins and assign the best scoring Human
					// Readable Description
					assignHumanReadableDescriptions();
					if (getSettings().hasGeneOntologyAnnotations() && getSettings().hasReferenceGoAnnotations()) {
						goAnnotsStringToObject();
					}
					// Evaluate AHRD's performance for each Protein:
					calculateEvaluationScores();
					// Estimate average performance of current Parameters:
					calcAveragesOfEvalScoreTPRandFPR();
//					if(getSettings().getParameters().getOrigin().equals("seed")) {
//						writeProteins(generation);
//					}
//				}
			}

			// Survival of the fittest
			NavigableSet<Parameters> fittnessRanking = new TreeSet<Parameters>();
			fittnessRanking.addAll(population);
			population.clear();
			while (fittnessRanking.size() > numberOfSurvivors) {
				fittnessRanking.pollFirst();
			}
			population.addAll(fittnessRanking);

			// Recombination of fit survivors
			Set<Set<Parameters>> uniqueMatingPairs = new HashSet<Set<Parameters>>();
			if (Utils.binomial(fittnessRanking.size(), 2) >= numberOfOffspring) { // Avoid endless loops occurring if size of fittnessRanking very small 
				while (population.size() < numberOfSurvivors + numberOfOffspring) {
					Set<Parameters> matingPair = new HashSet<Parameters>();
					do {
						matingPair.clear();
						while (matingPair.size() < 2) {
							matingPair.add(getRandomFitIndividual(fittnessRanking));
						}
					} while (uniqueMatingPairs.contains(matingPair));
					uniqueMatingPairs.add(matingPair);
					Parameters[] matingPairArray = matingPair.toArray(new Parameters[0]);
					population.add(matingPairArray[0].recombine(matingPairArray[1]));
				}
			}

			// Mutants of fit survivors
			while (population.size() < numberOfSurvivors + numberOfOffspring + numberOfMutants) {
				population.add(getRandomFitIndividual(fittnessRanking).neighbour(null));
			}

			// Fill the rest of the population with new parameter sets
			while (population.size() <= getSettings().getPopulationSize()) {
				population.add(Parameters.randomParameters(sortedDistinctBlastDatabaseNames));
			}

			// Remember the best parameter set and the generation it was found
			// in
			if (getBestParameters() != null) {
				diffAvgEvalScoreToLastGeneration = fittnessRanking.last().getAvgEvaluationScore() - getBestParameters().getAvgEvaluationScore();
			}
			if (getBestParameters() == null
					|| fittnessRanking.last().getAvgEvaluationScore() > getBestParameters().getAvgEvaluationScore()) {
				setBestParameters(fittnessRanking.last().clone());
				setGenerationBestParametersWereFoundIn(generation);
			}
			// Write output of current iteration:
			this.outWriter.writeIterationOutput(generation, getBestParameters(), diffAvgEvalScoreToLastGeneration, getBestParameters().getOrigin());
			generation += 1;
		}
	}
	
	/**
	 * Writes useful info for every protein according to the currents parameter set to file.
	 * Is meant for Debugging.
	 * Uses the generation number as file name.
	 * @param generation
	 * @throws IOException
	 */
	private void writeProteins(int generation) throws IOException{
		BufferedWriter outBufWrtr = new BufferedWriter(new FileWriter(generation + ".tsv"));
		outBufWrtr.write("QueryAccession\tBlastAccession\tSemSimGoAnnotationScore\tGOterms\tDescriptionScore\tLexicalScore\tRelativeBlastScore\tBitScore\tMaxBitScore\tDescription\n");
		int count = 0;
		for(Protein p:getProteins().values()) {
			if(p.getDescriptionScoreCalculator()==null) {
				System.out.println("NPE");
			} else {
				if (p.getDescriptionScoreCalculator().getHighestScoringBlastResult()==null) {
					//outBufWrtr.write(p.getAccession() + "\t-\n");
				} else {
					if(p.getDescriptionScoreCalculator().getHighestScoringBlastResult().getDescription()==null) {
						System.out.println("NPE");
					} else {
						count++;
						if (count == 178) {
							double test = p.getLexicalScoreCalculator().lexicalScore(p.getDescriptionScoreCalculator().getHighestScoringBlastResult());
							test = test + 0.0;
						}
						outBufWrtr.write(p.getAccession() + "\t"
					+ p.getDescriptionScoreCalculator().getHighestScoringBlastResult().getShortAccession() + "\t"
					+ p.getEvaluationScoreCalculator().getSemSimGoAnnotationScore() + "\t"
					+ Utils.joinStringCollection(",", p.getGoResults()) + "\t"
					+ p.getDescriptionScoreCalculator().getHighestScoringBlastResult().getDescriptionScore() + "\t"
					+ p.getLexicalScoreCalculator().lexicalScore(p.getDescriptionScoreCalculator().getHighestScoringBlastResult()) + "\t"
					+ p.getDescriptionScoreCalculator().relativeBlastScore(p.getDescriptionScoreCalculator().getHighestScoringBlastResult()) + "\t"
					+ p.getDescriptionScoreCalculator().getHighestScoringBlastResult().getBitScore() + "\t"
					+ p.getDescriptionScoreCalculator().getMaxBitScore() + "\t"
					+ p.getDescriptionScoreCalculator().getHighestScoringBlastResult().getDescription() + "\n");
					}
				}
			}
		}
		outBufWrtr.close();	
	}
	
	/**
	 * Returns a random individual (parameter set) from a navigable set of
	 * parameter sets ordered according to their fittness (evaluation score) A
	 * strong bias towards fitter individuals is applied
	 * 
	 * @param fittnessRanking
	 * @return random fit parameter set
	 */
	private static Parameters getRandomFitIndividual(NavigableSet<Parameters> fittnessRanking) {
		Parameters randomFitIndividual = null;
		int placeInRanking = Integer.MAX_VALUE;
		while (placeInRanking > fittnessRanking.size()) {
			placeInRanking = (int) Math.ceil(Math.abs(Utils.random.nextGaussian() * ((double) fittnessRanking.size() / 3)));
		}
		Iterator<Parameters> decendingRankingIter = fittnessRanking.descendingIterator();
		for (int i = 1; i <= placeInRanking; i++) {
			randomFitIndividual = decendingRankingIter.next();
		}
		return randomFitIndividual;
	}

	/**
	 * Calculates the average of AHRD's EvaluationScore (objective-function).
	 * If GO term scores have been computed the average is based upon them.
	 * Otherwise the conventional HRD based scores are used. If so also calculates 
	 * the average True-Positives- and False-Positives-Rates.
	 */
	public void calcAveragesOfEvalScoreTPRandFPR() {
		// average evaluation-score
		Double avgEvlScr = 0.0;
		// average TPR:
		Double avgTruePosRate = 0.0;
		// average FPR:
		Double avgFalsePosRate = 0.0;
		// Evaluate GO annotations.
		if (getSettings().hasGeneOntologyAnnotations() && getSettings().hasReferenceGoAnnotations()) {
			for (Protein p : getProteins().values()) {
				EvaluationScoreCalculator e = p.getEvaluationScoreCalculator();
				if (e != null) {
					//Depending on the settings the go annotation f-score with the highest level of complexity is used
					if (getSettings().getCalculateSemSimGoF1Scores()) {
						avgEvlScr += e.getSemSimGoAnnotationScore();
					} else {
						if (getSettings().getCalculateAncestryGoF1Scores()) {
							avgEvlScr += e.getAncestryGoAnnotationScore();
						} else {
							avgEvlScr += e.getSimpleGoAnnotationScore();
						}
					}
				}
			}
		} else { // Otherwise use HRD based scores
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
		}
		// average each number:
		Double numberOfProts = new Double(getProteins().size());
		avgEvlScr = avgEvlScr / numberOfProts;
		avgTruePosRate = avgTruePosRate / numberOfProts;
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
