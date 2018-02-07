package ahrd.controller;

import static ahrd.controller.Settings.getSettings;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import ahrd.exception.MissingAccessionException;
import ahrd.exception.MissingInterproResultException;
import ahrd.model.BlastResult;
import ahrd.model.DescriptionScoreCalculator;
import ahrd.model.EvaluationScoreCalculator;
import ahrd.model.Fscore;
import ahrd.model.GOterm;
import ahrd.model.Protein;
import ahrd.model.TokenScoreCalculator;
import ahrd.view.TrainerOutputWriter;

public abstract class Trainer extends Evaluator {

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
	 * Uses the currently assigned set of parameters to annotate the proteins with either GO terms or Human Readable Descriptions.
	 * From the annotation an Fscore as well as Precision and Recall are derived by comparison to a ground truth.
	 * If term centric annotation is requested is will be derived from the protein centric go annotation and the scores are calculated based on the term centric annotation
	 * @throws MissingInterproResultException
	 * @throws IOException
	 * @throws SQLException
	 * @throws MissingAccessionException
	 */
	protected void scoreCurrentParameters() throws MissingInterproResultException, IOException, SQLException, MissingAccessionException {
		reinitializeBlastResults();
		if (getSettings().hasGeneOntologyAnnotations() && getSettings().hasGroundTruthGoAnnotations()) {
			assignGeneOntologyTerms(); // Iterate over Proteins and assign go terms
			if(getSettings().hasGoTermCentricTermsFile()){
				termCentricAnnotation(); // Use Protein centric go annotations to derive term centric annotations (with confidence coefficients)
			} else {
				goAnnotsStringToObject(); // Needed only when evaluation is to be based on protein centric go annotations  
				calculateEvaluationScores(); // Evaluate AHRD's performance for each Protein
			}
		} else {
			assignHumanReadableDescriptions(); // Iterate over all Proteins and assign the best scoring Human Readable Description
			calculateEvaluationScores(); // Evaluate AHRD's performance for each Protein
		}
		// Estimate average performance of current Parameters:
		calcAveragesOfEvalScorePrecisionAndRecall();
	}

	/**
	 * Calculates the average of AHRD's EvaluationScore (objective-function).
	 * If GO term scores have been computed the average is based upon them.
	 * Further, if the GO annotations have been requested to be performed in a term centric manner the average scores are based upon them.  
	 * Otherwise the conventional HRD based scores are used.
	 */
	public void calcAveragesOfEvalScorePrecisionAndRecall() {
		Double numberOfProts = new Double(getProteins().size());
		// average evaluation-score
		Double avgEvlScr = 0.0;
		// average Precision (PPV):
		Double avgPrecision = 0.0;
		// average Recall (TPR):
		Double avgRecall = 0.0;
		// Evaluate GO annotations.
		if (getSettings().hasGeneOntologyAnnotations() && getSettings().hasGroundTruthGoAnnotations()) {
			if (getSettings().hasGoTermCentricTermsFile()) { // Evaluate term centric go annotations
				Map<String, Fscore> bestTermCentricScores = new HashMap<String, Fscore>();
				for (GOterm term : this.goCentricTerms) {
					bestTermCentricScores.put(term.getAccession(), new Fscore());
				}
				// Iterate over annotation confidence thresholds
				for(int i=0; i<=100; i++) {
					double threshold = ((double) i) * 0.01;
					Map<String, Fscore> scores = this.calcTermCentricScores(threshold);
					for (String termAcc : scores.keySet()) {
						if (scores.get(termAcc).getScore() > bestTermCentricScores.get(termAcc).getScore()) {
							bestTermCentricScores.put(termAcc, scores.get(termAcc));
						}
					}
				}
				// Average over all term centric terms
				for (String termAcc : bestTermCentricScores.keySet()) {
					avgEvlScr += bestTermCentricScores.get(termAcc).getScore();
					avgPrecision += bestTermCentricScores.get(termAcc).getPrecision();
					avgRecall += bestTermCentricScores.get(termAcc).getRecall();
				}
				avgEvlScr = avgEvlScr / (double) bestTermCentricScores.size();
				avgPrecision = avgPrecision / (double) bestTermCentricScores.size();
				avgRecall = avgRecall / (double) bestTermCentricScores.size();
			} else { // Evaluate protein centric go annotations 
				for (Protein p : getProteins().values()) {
					EvaluationScoreCalculator e = p.getEvaluationScoreCalculator();
					if (e != null) {
						//Depending on the settings the go annotation f-score with the highest level of complexity is used
						if (getSettings().doCalculateSemSimGoF1Scores()) {
							avgEvlScr += e.getSemSimGoAnnotationScore().getScore();
							avgPrecision += e.getSemSimGoAnnotationScore().getPrecision();
							avgRecall += e.getSemSimGoAnnotationScore().getRecall();
						} else {
							if (getSettings().doCalculateAncestryGoF1Scores()) {
								avgEvlScr += e.getAncestryGoAnnotationScore().getScore();
								avgPrecision += e.getAncestryGoAnnotationScore().getPrecision();
								avgRecall += e.getAncestryGoAnnotationScore().getRecall();
							} else {
								// getSettings().doCalculateSimpleGoF1Scores() Needs to be checked?
								avgEvlScr += e.getSimpleGoAnnotationScore().getScore();
								avgPrecision += e.getSimpleGoAnnotationScore().getPrecision();
								avgRecall += e.getSimpleGoAnnotationScore().getRecall();
							}
						}
					}
				}
				// average each number:
				avgEvlScr = avgEvlScr / numberOfProts;
				avgPrecision = avgPrecision / numberOfProts;
				avgRecall = avgRecall / numberOfProts;
			}
		} else { // Otherwise use HRD based scores
			for (Protein p : getProteins().values()) {
				EvaluationScoreCalculator e = p.getEvaluationScoreCalculator();
				if (e != null) {
					if (e.getEvalutionScore() != null) {
						avgEvlScr += e.getEvalutionScore().getScore();
						avgPrecision += e.getEvalutionScore().getPrecision();
						avgRecall += e.getEvalutionScore().getRecall();
					}
				}
			}
			// average each number:
			avgEvlScr = avgEvlScr / numberOfProts;
			avgPrecision = avgPrecision / numberOfProts;
			avgRecall = avgRecall / numberOfProts;
		}
		// Annotate the current settings with the averaged scores:
		getSettings().setAvgEvaluationScore(avgEvlScr);
		getSettings().setAvgPrecision(avgPrecision);
		getSettings().setAvgRecall(avgRecall);
	}

	/**
	 * This calculates the average maximum evaluation score AHRD could possibly
	 * achieve.
	 */
	public void calcAvgMaxEvaluationScore() {
		double avgMaxEvlScr = 0.0; 		// init average maximum evaluation-score
		if (getSettings().hasGeneOntologyAnnotations() && getSettings().hasGroundTruthGoAnnotations()) { 		// Evaluate GO annotations.
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
				prot.getEvaluationScoreCalculator().findBlastResultWithHighestPossibleDescriptionScore();
				avgMaxEvlScr += prot.getEvaluationScoreCalculator().getHighestPossibleDescriptionScore().getScore();
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
