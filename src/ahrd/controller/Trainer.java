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

public class Trainer extends Evaluator {

	protected Parameters bestParameters;
	protected TrainerOutputWriter outWriter;

	/**
	 * The average of AHRD's maximum evaluation score for each Protein. This is
	 * the maximum of the evaluation scores calculated for all Descriptions of
	 * each Protein. These maximums are then averaged.
	 */
	protected Double avgMaxEvaluationScore = 0.0;

	/**
	 * Constructor initializes the Settings as given in the argument input.yml
	 * 
	 * @param pathToInputYml
	 * @throws IOException
	 */
	public Trainer(String pathToInputYml) throws IOException {
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
	 * Calculates the average of AHRD's EvaluationScore (objective-function).
	 * If GO term scores have been computed the average is based upon them.
	 * Otherwise the conventional HRD based scores are used.
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
					if (getSettings().doCalculateSemSimGoF1Scores()) {
						avgEvlScr += e.getSemSimGoAnnotationScore().getScore();
					} else {
						if (getSettings().doCalculateAncestryGoF1Scores()) {
							avgEvlScr += e.getAncestryGoAnnotationScore().getScore();
						} else {
							avgEvlScr += e.getSimpleGoAnnotationScore().getScore();
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
		double avgMaxEvlScr = 0.0; 		// init average maximum evaluation-score
		if (getSettings().hasGeneOntologyAnnotations() && getSettings().hasReferenceGoAnnotations()) { 		// Evaluate GO annotations.
			for (Protein p : getProteins().values()) {
				EvaluationScoreCalculator e = p.getEvaluationScoreCalculator();
				e.findHighestPossibleGoScore();
				//Depending on the settings the go annotation f-score with the highest level of complexity is used
				if (getSettings().doCalculateSemSimGoF1Scores()) {
					avgMaxEvlScr += e.getHighestPossibleSemSimGoAnnotationScore().getScore();
				} else {
					if (getSettings().doCalculateAncestryGoF1Scores()) {
						avgMaxEvlScr += e.getHighestPossibleAncestryGoAnnotationScore().getScore();
					} else {
						avgMaxEvlScr += e.getHighestPossibleSimpleGoAnnotationScore().getScore();
					}
				}
			}
		} else { // Otherwise use HRD based scores
			for (Protein prot : getProteins().values()) {
				prot.getEvaluationScoreCalculator().findHighestPossibleEvaluationScore();
				avgMaxEvlScr += prot.getEvaluationScoreCalculator().getHighestPossibleEvaluationScore();
			}
		}
		setAvgMaxEvaluationScore(avgMaxEvlScr / getProteins().size());		// calculate average
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
		+ getSettings().getAvgEvaluationScore() + "\t"
		+ getSettings().getParameters().getOrigin() + "\t" 
		+ getSettings().getTokenScoreBitScoreWeight() + "\t"
		+ getSettings().getTokenScoreDatabaseScoreWeight() + "\t"
		+ getSettings().getTokenScoreOverlapScoreWeight());
		for (String blastDb : sortedBlastDatabases) {
			outBufWrtr.write("\t" + getSettings().getBlastDbWeight(blastDb));
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

	public void setBestParameters(Parameters bestParameters) {
		this.bestParameters = bestParameters;
	}

	public Double getAvgMaxEvaluationScore() {
		return avgMaxEvaluationScore;
	}

	public void setAvgMaxEvaluationScore(Double avgMaxEvaluationScore) {
		this.avgMaxEvaluationScore = avgMaxEvaluationScore;
	}

}
