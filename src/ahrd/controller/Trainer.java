package ahrd.controller;

import static ahrd.controller.Settings.getSettings;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import ahrd.model.BlastResult;
import ahrd.model.DescriptionScoreCalculator;
import ahrd.model.EvaluationScoreCalculator;
import ahrd.model.Protein;
import ahrd.model.TokenScoreCalculator;
import ahrd.view.TrainerOutputWriter;

public abstract class Trainer extends Evaluator {

	protected Parameters bestParameters;
	protected TrainerOutputWriter outputWriter;

	/**
	 * The average of AHRD's maximum description evaluation score for each protein.
	 * This is the maximum of the evaluation scores calculated for all descriptions of each Protein.
	 * These maximums are then averaged.
	 */
	protected Double avgMaxDescriptionScore = 0.0;
	/**
	 * The average of AHRD's maximum GO evaluation score for each protein.
	 * This is the maximum of the evaluation scores calculated for all GO annotation of each Protein.
	 * These maximums are then averaged.
	 */
	protected Double avgMaxGoScore = 0.0;

	/**
	 * Constructor initializes the Settings as given in the argument input.yml
	 * 
	 * @param pathToInputYml
	 * @throws IOException
	 */
	protected Trainer(String pathToInputYml) throws IOException {
		super(pathToInputYml);
	}

	/**
	 * Removes all BlastResults from each Protein and adds them again.
	 * Reinitializes the TokenScoreCalculator and the DescriptionScoreCalculator of each Protein.
	 * This triggers recalculation of cumulative and total token scores in the TokenScoreCalulator and the maxBitScore in the DescriptionScoreCalculator.
	 * The target of this operation are the CumulativeTokenBlastDatabaseScore and the TotalTokenBlastDatabaseScore in the TokenScoreCalculator.
	 * These scores depend on the BlastDbWeights and can't be transfered from other sets of Parameters.  
	 */
	protected void reinitializeBlastResults() {
		for (Protein p : getProteins().values()) {
			Map<String, List<BlastResult>> blastResults = p.getBlastResults();
			p.setBlastResults(new HashMap<String, List<BlastResult>>());
			p.setTokenScoreCalculator(new TokenScoreCalculator(p));
			p.setDescriptionScoreCalculator(new DescriptionScoreCalculator(p));
			for(String blastDb : blastResults.keySet()) {
				for(BlastResult br : blastResults.get(blastDb)){
					p.addBlastResult(br);
				}
			}
		}
	}

	/**
	 * Calculates the average of AHRD's description evaluation score.
	 * In the description parameter optimization (training) this is used as objective-function.
	 * The average recall and precision are calculated as well 
	 */
	public void calcAveragesOfDescriptionScores() {
		// average evaluation-score
		Double avgEvlScr = 0.0;
		// average Precision (PPV):
		Double avgPrecision = 0.0;
		// average Recall (TPR):
		Double avgRecall = 0.0;
		// Number of proteins covert by prediction
		Integer numberOfProts = 0;
		for (Protein p : getProteins().values()) {
			EvaluationScoreCalculator e = p.getEvaluationScoreCalculator();
			if (e != null && e.getEvalutionScore() != null && !e.getEvalutionScore().getScore().isNaN()) {
					avgEvlScr += e.getEvalutionScore().getScore();
					avgPrecision += e.getEvalutionScore().getPrecision();
					avgRecall += e.getEvalutionScore().getRecall();
					numberOfProts++;
			}
		}
		// average each number:
		if (avgEvlScr > 0.0)
			avgEvlScr = avgEvlScr / numberOfProts;
		if (avgPrecision > 0.0)
			avgPrecision = avgPrecision / numberOfProts;
		if (avgRecall > 0.0)
			avgRecall = avgRecall / numberOfProts;
		// done:
		getSettings().setAvgDescriptionEvaluationScore(avgEvlScr);
		getSettings().setAvgDescriptionPrecision(avgPrecision);
		getSettings().setAvgDescriptionRecall(avgRecall);
	}

	/**
	 * Calculates the average of AHRD's GO evaluation score.
	 * AHRD can calculate three different GO scores.
	 * Which GO scores are calculated tepends on AHRD's settings.
	 * If multiple GO scores are requested to be calculated the one with the highest complexity is used here.
	 * In the GO parameter optimization (training) this is used as objective-function.
	 * The average recall and precision are calculated as well.
	 */
	public void calcAveragesOfGoScores() {
		// average evaluation-score
		Double avgEvlScr = 0.0;
		// average Precision (PPV):
		Double avgPrecision = 0.0;
		// average Recall (TPR):
		Double avgRecall = 0.0;
		// Number of proteins covert by prediction
		Integer numberOfProts = 0;
		// Evaluate GO annotations.
		for (Protein p : getProteins().values()) {
			EvaluationScoreCalculator e = p.getEvaluationScoreCalculator();
			if (e != null) {
				//Depending on the settings the go annotation f-score with the highest level of complexity is used
				if (getSettings().doCalculateSemSimGoF1Scores()) {
					if (!e.getSemSimGoAnnotationScore().getScore().isNaN()) {
						avgEvlScr += e.getSemSimGoAnnotationScore().getScore();
						avgPrecision += e.getSemSimGoAnnotationScore().getPrecision();
						avgRecall += e.getSemSimGoAnnotationScore().getRecall();
						numberOfProts++;							
					}
				} else {
					if (getSettings().doCalculateAncestryGoF1Scores()) {
						if (!e.getAncestryGoAnnotationScore().getScore().isNaN()) {
							avgEvlScr += e.getAncestryGoAnnotationScore().getScore();
							avgPrecision += e.getAncestryGoAnnotationScore().getPrecision();
							avgRecall += e.getAncestryGoAnnotationScore().getRecall();
							numberOfProts++;
						}
					} else {
						if (!e.getSimpleGoAnnotationScore().getScore().isNaN()) {
							avgEvlScr += e.getSimpleGoAnnotationScore().getScore();
							avgPrecision += e.getSimpleGoAnnotationScore().getPrecision();
							avgRecall += e.getSimpleGoAnnotationScore().getRecall();
							numberOfProts++;
						}
					}
				}
			}
		}
		// average each number:
		if (avgEvlScr > 0.0)
			avgEvlScr = avgEvlScr / numberOfProts;
		if (avgPrecision > 0.0)
			avgPrecision = avgPrecision / numberOfProts;
		if (avgRecall > 0.0)
			avgRecall = avgRecall / numberOfProts;
		// done:
		getSettings().setAvgGoEvaluationScore(avgEvlScr);
		getSettings().setAvgGoPrecision(avgPrecision);
		getSettings().setAvgGoRecall(avgRecall);
	}

	
	/**
	 * This calculates the average maximum description evaluation score AHRD could possibly achieve.
	 */
	public void calcAvgMaxDescriptionScore() {
		double avgMaxEvlScr = 0.0; 		// init average maximum evaluation-score
		Integer numberOfProts = 0;
		for (Protein prot : getProteins().values()) {
			prot.getEvaluationScoreCalculator().findBlastResultWithHighestPossibleDescriptionScore();
			if (!prot.getEvaluationScoreCalculator().getHighestPossibleDescriptionScore().getScore().isNaN()) {
				avgMaxEvlScr += prot.getEvaluationScoreCalculator().getHighestPossibleDescriptionScore().getScore();
				numberOfProts++;
			}
		}
		setAvgMaxDescriptionScore(avgMaxEvlScr / numberOfProts);		// calculate average
	}
		
	
	/**
	 * This calculates the average maximum GO evaluation score AHRD could possibly achieve.
	 */
	public void calcAvgMaxGoScore() {
		double avgMaxEvlScr = 0.0; 		// init average maximum evaluation-score
		Integer numberOfProts = 0;
		for (Protein p : getProteins().values()) {
			EvaluationScoreCalculator e = p.getEvaluationScoreCalculator();
			e.findHighestPossibleGoScore();
			//Depending on the settings the go annotation f-score with the highest level of complexity is used
			if (getSettings().doCalculateSemSimGoF1Scores()) {
				if (!e.getHighestPossibleSemSimGoAnnotationScore().getScore().isNaN()) {
					avgMaxEvlScr += e.getHighestPossibleSemSimGoAnnotationScore().getScore();
					numberOfProts++;
				}
			} else {
				if (getSettings().doCalculateAncestryGoF1Scores()) {
					if (!e.getHighestPossibleAncestryGoAnnotationScore().getScore().isNaN()) {
						avgMaxEvlScr += e.getHighestPossibleAncestryGoAnnotationScore().getScore();
						numberOfProts++;
					}
				} else {
					if (!e.getHighestPossibleSimpleGoAnnotationScore().getScore().isNaN()) {
						avgMaxEvlScr += e.getHighestPossibleSimpleGoAnnotationScore().getScore();
						numberOfProts++;
					}
				}
			}
		}
		setAvgMaxGoScore(avgMaxEvlScr / numberOfProts);		// calculate average
	}

	/**
	 * Writes useful info for every protein according to the currents parameter set to file.
	 * Is meant for Debugging.
	 * Uses the generation number as file name.
	 * @param iteration
	 * @throws IOException
	 */
	protected void writeProteins(int iteration) throws IOException{
		BufferedWriter outBufWrtr = new BufferedWriter(new FileWriter(iteration + ".tsv"));
		outBufWrtr.write("#Generation\tAverage Evaluation-Score(F-Score)\tOrigin\tToken-Score-Bit-Score-Weight\tToken-Score-Database-Score-Weight\tToken-Score-Overlap-Score-Weight");
		List<String> sortedBlastDatabases = new ArrayList<String>(getSettings().getBlastDatabases());
		Collections.sort(sortedBlastDatabases);
		for (String blastDb : sortedBlastDatabases) {
			outBufWrtr.write("\t" + blastDb + "-Weight");
			outBufWrtr.write("\t" + blastDb + "-Description-Score-Bit-Score-Weight");
		}
		outBufWrtr.write("\n");
		outBufWrtr.write("#" + iteration + "\t" 
		+ getSettings().getAvgDescriptionEvaluationScore() + "\t"
//		+ getSettings().getParameters().getOrigin() + "\t" 
		+ getSettings().getDescriptionTokenScoreBitScoreWeight() + "\t"
		+ getSettings().getDescriptionTokenScoreDatabaseScoreWeight() + "\t"
		+ getSettings().getDescriptionTokenScoreOverlapScoreWeight());
		for (String blastDb : sortedBlastDatabases) {
			outBufWrtr.write("\t" + getSettings().getDescriptionBlastDbWeight(blastDb));
			outBufWrtr.write("\t" + getSettings().getDescriptionScoreBitScoreWeight(blastDb));
		}
		outBufWrtr.write("\n");
		
		outBufWrtr.write("QueryAccession\tBlastAccession\tSemSimGoAnnotationScore\tGOterms\tDescriptionScore\tLexicalScore\tRelativeBlastScore\tBitScore\tMaxBitScore\tDescription\tTokenScoreHighScore\tTotalTokenBitScore\tTotalTokenBlastDatabaseScore\tTotalTokenOverlapScore\tCumulativeTokenBitScores\tCumulativeTokenBlastDatabaseScores\tgetCumulativeTokenOverlapScores\n");
		for(Protein p:getProteins().values()) {
			if (p.getDescriptionScoreCalculator().getHighestScoringBlastResult() == null) {
				outBufWrtr.write(p.getAccession() + "\t\t" 
			+ p.getEvaluationScoreCalculator().getSemSimGoAnnotationScore() + "\t"
			+ Utils.joinStringCollection(",", p.getGoResults()) + "\t"
			+ "0\t0\t0\t0\t0\t\t0\t0\n");
			} else {
				if(p.getDescriptionScoreCalculator().getHighestScoringBlastResult().getDescription() == null) {
					System.out.println("NPE");
				} else {
					outBufWrtr.write(p.getAccession() + "\t"
				+ p.getDescriptionScoreCalculator().getHighestScoringBlastResult().getShortAccession() + "\t"
				+ p.getEvaluationScoreCalculator().getSemSimGoAnnotationScore() + "\t"
				+ Utils.joinStringCollection(",", p.getGoResults()) + "\t"
				+ p.getDescriptionScoreCalculator().getHighestScoringBlastResult().getDescriptionScore() + "\t"
				+ p.getLexicalScoreCalculator().lexicalScore(p.getDescriptionScoreCalculator().getHighestScoringBlastResult()) + "\t"
				+ p.getDescriptionScoreCalculator().relativeBlastScore(p.getDescriptionScoreCalculator().getHighestScoringBlastResult()) + "\t"
				+ p.getDescriptionScoreCalculator().getHighestScoringBlastResult().getBitScore() + "\t"
				+ p.getDescriptionScoreCalculator().getMaxBitScore() + "\t"
				+ p.getDescriptionScoreCalculator().getHighestScoringBlastResult().getDescription() + "\t"
				+ p.getTokenScoreCalculator().getTokenHighScore() + "\t"
				+ p.getTokenScoreCalculator().getTotalTokenBitScore() + "\t"
				+ p.getTokenScoreCalculator().getTotalTokenBlastDatabaseScore() + "\t"
				+ p.getTokenScoreCalculator().getTotalTokenOverlapScore() + "\t"
				+ p.getTokenScoreCalculator().getCumulativeTokenBitScores() + "\t"
				+ p.getTokenScoreCalculator().getCumulativeTokenBlastDatabaseScores() + "\t"
				+ p.getTokenScoreCalculator().getCumulativeTokenOverlapScores()
				+ "\n");
				}
			}
		}
		outBufWrtr.close();
		
		File theDir = new File(Integer.toString(iteration));
		theDir.mkdir();
		for(Protein p:getProteins().values()) {
			outBufWrtr = new BufferedWriter(new FileWriter("./" + iteration + "/" + p.getAccession() + ".tsv"));
			outBufWrtr.write("BlastDb\tBlastAccession\tShortAccession\tDescriptionScore\tDescription\tBitScore\tTokens\tInformativeTokens\tTokenScores\tDescriptionScoreRelativeBlastScore\tLexicalScore\tLexicalScoreCorrectionFactor\tLexicalScoreSumTokenScoresDividebByHighScore\tSumOfAllTokenScores\n");
			for(String blastDb : p.getBlastResults().keySet()) {
				for(BlastResult br : p.getBlastResults().get(blastDb)){
					TokenScoreCalculator tsk = p.getTokenScoreCalculator();
					Set<String> informativeTokens = new HashSet<String>();
					List<String> tokenScores = new ArrayList<String>();
					for(String token : br.getTokens()) {
						if (tsk.isInformativeToken(token)) {
							informativeTokens.add(token);
						}
						tokenScores.add(tsk.getTokenScores().get(token).toString());
					}
					outBufWrtr.write(blastDb + "\t"
							+ br.getAccession() + "\t"
							+ br.getShortAccession() + "\t" 
							+ br.getDescriptionScore() + "\t"
							+ br.getDescription() + "\t"
							+ br.getBitScore() + "\t"
							+ Utils.joinStringCollection(",", br.getTokens()) + "\t"
							+ Utils.joinStringCollection(",", informativeTokens) + "\t"
							+ Utils.joinStringCollection(",", tokenScores) + "\t"
							+ p.getDescriptionScoreCalculator().relativeBlastScore(br) + "\t"
							+ p.getLexicalScoreCalculator().lexicalScore(br) + "\t"
							+ p.getLexicalScoreCalculator().correctionFactor(br) + "\t"
							+ p.getLexicalScoreCalculator().sumTokenScoresDividebByHighScore(br) + "\t"
							+ p.getTokenScoreCalculator().sumOfAllTokenScores(br) + "\t"
							+ "\n");
				}
			}
			outBufWrtr.close();
		}
	}
	
	public Parameters getBestParameters() {
		return bestParameters;
	}

	public void setBestParameters(Parameters parameters) {
		this.bestParameters = parameters;
	}
	
	public Double getAvgMaxDescriptionScore() {
		return avgMaxDescriptionScore;
	}

	public void setAvgMaxDescriptionScore(Double score) {
		this.avgMaxDescriptionScore = score;
	}
	
	public Double getAvgMaxGoScore() {
		return avgMaxGoScore;
	}

	public void setAvgMaxGoScore(Double score) {
		this.avgMaxGoScore = score;
	}

	public TrainerOutputWriter getOutputWriter() {
		return outputWriter;
	}

	public void setOutputWriter(TrainerOutputWriter outputWriter) {
		this.outputWriter = outputWriter;
	}

}
