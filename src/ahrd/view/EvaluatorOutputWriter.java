package ahrd.view;

import static ahrd.controller.Settings.getSettings;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;

import ahrd.controller.AHRD;
import ahrd.model.BlastResult;
import ahrd.model.CompetitorAnnotation;
import ahrd.model.Fscore;
import ahrd.model.Protein;

public class EvaluatorOutputWriter extends TsvOutputWriter {

	protected BufferedWriter hrdScoresWriter;

	public EvaluatorOutputWriter(Collection<Protein> proteins) {
		super(proteins);
	}

	public void writeOutput() throws IOException {
		BufferedWriter bw = new BufferedWriter(new FileWriter(getSettings().getPathToOutput()));
		if (getSettings().doWriteHRDScoresToOutput())
			writeHRDScoresOutputHeader();

		// Header
		bw.write("# AHRD-Version " + AHRD.VERSION + " (Evaluator)\n");
		bw.write("\n");
		// Column-Names:
		bw.write(ahrdColumnNames());
		bw.write("\tHRD-Length\tGround-Truth-Description\tRef-Lenght\tEvaluation-Score");
		if (getSettings().doWriteFscoreDetailsToOutput()) {
			bw.write("\tDiff-to-bestCompetitor\tPrecision\tRecall\tRemaining-Uncertainty\tMissing-Information");
		}
		if (getSettings().getWriteTokenSetToOutput()) {
			bw.write("\t\"Tokens (tkn->score)\"");
		}
		if (getSettings().getWriteScoresToOutput()) {
			bw.write("\tSum(Token-Scores)\tTokenHighScore\tCorrection-Factor\tLexical-Score\tRelativeBitScore");
		}
		if (getSettings().getWriteBestBlastHitsToOutput()) {
			bw.write(buildBestBlastHitsHeader());
		}
		if (getSettings().hasCompetitors()) {
			for (String competitor : getSettings().getCompetitorSettings().keySet()) {
				bw.write("\t" + competitor + "-Description\t" + competitor + "-Description-Length\t" + competitor + "-Description-Evaluation-Score");
				if (getSettings().doWriteFscoreDetailsToOutput()) {
					bw.write("\t" + competitor + "-Description-Evaluation-Score-Precision");
					bw.write("\t" + competitor + "-Description-Evaluation-Score-Recall");
					bw.write("\t" + competitor + "-Description-Evaluation-Score-Remaining-Uncertainty");
					bw.write("\t" + competitor + "-Description-Evaluation-Score-Missing-Information");
				}
				if (getSettings().hasGeneOntologyAnnotations() && getSettings().hasGroundTruthGoAnnotations()) {
					bw.write("\t" + competitor + "-GO-Annotations");
					if (getSettings().doCalculateSimpleGoF1Scores()){
						bw.write("\t" + competitor + "-GO-Annotations-Simple-F-Score");
						if (getSettings().doWriteFscoreDetailsToOutput()) {
							bw.write("\t" + competitor + "-GO-Annotations-Simple-F-Score-Precision");
							bw.write("\t" + competitor + "-GO-Annotations-Simple-F-Score-Recall");
							bw.write("\t" + competitor + "-GO-Annotations-Simple-F-Score-Remaining-Uncertainty");
							bw.write("\t" + competitor + "-GO-Annotations-Simple-F-Score-Missing-Information");
						}
					}
					if (getSettings().doCalculateAncestryGoF1Scores()){
						bw.write("\t" + competitor + "-GO-Annotations-Ancestry-F-Score");
						if (getSettings().doWriteFscoreDetailsToOutput()) {
							bw.write("\t" + competitor + "-GO-Annotations-Ancestry-F-Score-Precision");
							bw.write("\t" + competitor + "-GO-Annotations-Ancestry-F-Score-Recall");
							bw.write("\t" + competitor + "-GO-Annotations-Ancestry-F-Score-Remaining-Uncertainty");
							bw.write("\t" + competitor + "-GO-Annotations-Ancestry-F-Score-Missing-Information");
						}
					}
					if (getSettings().doCalculateSemSimGoF1Scores()){
						bw.write("\t" + competitor + "-GO-Annotations-SemSim-F-Score");
						if (getSettings().doWriteFscoreDetailsToOutput()) {
							bw.write("\t" + competitor + "-GO-Annotations-SemSim-F-Score-Precision");
							bw.write("\t" + competitor + "-GO-Annotations-SemSim-F-Score-Recall");
							bw.write("\t" + competitor + "-GO-Annotations-SemSim-F-Score-Remaining-Uncertainty");
							bw.write("\t" + competitor + "-GO-Annotations-SemSim-F-Score-Missing-Information");						}
					}
				}
			}
		}
		if (getSettings().doFindHighestPossibleEvaluationScore()) {
			bw.write("\tBlast-Result-With-Highest-Possible-Description-Score");
			bw.write("\tHighest-Possible-Description-Score");
			if (getSettings().doWriteFscoreDetailsToOutput()) {
				bw.write("\tHighest-Possible-Description-Score-Precision");
				bw.write("\tHighest-Possible-Description-Score-Recall");
				bw.write("\tHighest-Possible-Description-Score-Remaining-Uncertainty");
				bw.write("\tHighest-Possible-Description-Score-Missing-Information");
			}
		}
		if (getSettings().hasGeneOntologyAnnotations() && getSettings().hasGroundTruthGoAnnotations()) {
			bw.write("\tGround-Truth-GO-Annotations");
			if (getSettings().doCalculateSimpleGoF1Scores()) {
				bw.write("\tAHRD-GO-Annotations-Simple-F-Score");
				if (getSettings().doWriteFscoreDetailsToOutput()) {
					bw.write("\tAHRD-GO-Annotations-Simple-F-Score-Precision\tAHRD-GO-Annotations-Simple-F-Score-Recall\tAHRD-GO-Annotations-Simple-F-Score-Remaining-Uncertainty\tAHRD-GO-Annotations-Simple-F-Score-Missing-Information");
				}
			}
			if (getSettings().doCalculateAncestryGoF1Scores()) {
				bw.write("\tAHRD-GO-Annotations-Ancestry-F-Score");
				if (getSettings().doWriteFscoreDetailsToOutput()) {
					bw.write("\tAHRD-GO-Annotations-Ancestry-F-Score-Precision\tAHRD-GO-Annotations-Ancestry-F-Score-Recall\tAHRD-GO-Annotations-Ancestry-F-Score-Remaining-Uncertainty\tAHRD-GO-Annotations-Ancestry-F-Score-Missing-Information");
				}
			}
			if (getSettings().doCalculateSemSimGoF1Scores()) {
				bw.write("\tAHRD-GO-Annotations-SemSim-F-Score");
				if (getSettings().doWriteFscoreDetailsToOutput()) {
					bw.write("\tAHRD-GO-Annotations-SemSim-F-Score-Precision\tAHRD-GO-Annotations-SemSim-F-Score-Recall\tAHRD-GO-Annotations-SemSim-F-Score-Remaining-Uncertainty\tAHRD-GO-Annotations-SemSim-F-Score-Missing-Information");
				}
			}
			if (getSettings().doFindHighestPossibleGoScore()) {
				if (getSettings().doCalculateSimpleGoF1Scores()) {
					bw.write("\tHighest-Possible-BlastResult-GO-Annotations-Simple-F-Score");
					if (getSettings().doWriteFscoreDetailsToOutput()) {
						bw.write("\tHighest-Possible-BlastResult-GO-Annotations-Simple-F-Score-Precision\tHighest-Possible-BlastResult-GO-Annotations-Simple-F-Score-Recall\tHighest-Possible-BlastResult-GO-Annotations-Simple-F-Score-Remaining-Uncertainty\tHighest-Possible-BlastResult-GO-Annotations-Simple-F-Score-Missing-Information");
					}
				}
				if (getSettings().doCalculateAncestryGoF1Scores()) {
					bw.write("\tHighest-Possible-BlastResult-GO-Annotations-Ancestry-F-Score");
					if (getSettings().doWriteFscoreDetailsToOutput()) {
						bw.write("\tHighest-Possible-BlastResult-GO-Annotations-Ancestry-F-Score-Precision\tHighest-Possible-BlastResult-GO-Annotations-Ancestry-F-Score-Recall\tHighest-Possible-BlastResult-GO-Annotations-Ancestry-F-Score-Remaining-Uncertainty\tHighest-Possible-BlastResult-GO-Annotations-Ancestry-F-Score-Missing-Information");
					}
				}
				if (getSettings().doCalculateSemSimGoF1Scores()) {
					bw.write("\tHighest-Possible-BlastResult-GO-Annotations-SemSim-F-Score");
					if (getSettings().doWriteFscoreDetailsToOutput()) {
						bw.write("\tHighest-Possible-BlastResult-GO-Annotations-SemSim-F-Score-Precision\tHighest-Possible-BlastResult-GO-Annotations-SemSim-F-Score-Recall\tHighest-Possible-BlastResult-GO-Annotations-SemSim-F-Score-Remaining-Uncertainty\tHighest-Possible-BlastResult-GO-Annotations-SemSim-F-Score-Missing-Information");
					}
				}
			}
		}
		bw.write("\n");

		for (Protein prot : getProteins()) {
			// Generate the Human Readable Description:
			String csvRow = buildDescriptionLine(prot, "\t");
			// Write out the Evaluator-Score and the Ground-Truth-Description:
			csvRow += buildDescriptionEvaluationColumns(prot);
			// Append further information, if requested:
			if (getSettings().getWriteTokenSetToOutput()) {
				csvRow += buildTokenSetCell(prot);
			}
			if (getSettings().getWriteScoresToOutput()) {
				csvRow += buildDescScoreCells(prot);
			}
			if (getSettings().getWriteBestBlastHitsToOutput()) {
				csvRow += buildBestBlastHitsColumns(prot);
			}
			if (getSettings().hasCompetitors()) {
				for (String competitor : getSettings().getCompetitorSettings().keySet()) {
					csvRow += buildCompetitorColumns(competitor, prot);
				}
			}
			if (getSettings().doFindHighestPossibleEvaluationScore()) {
				csvRow += buildBlastResultWithHighestPossibleDescriptionScoreColumns(prot);
			}
			if (getSettings().hasGeneOntologyAnnotations() && getSettings().hasGroundTruthGoAnnotations()) {
				csvRow += buildGroundTruthGoAnnotationColumns(prot);
			}
			if (getSettings().doFindHighestPossibleGoScore()) {
				if (getSettings().doCalculateSimpleGoF1Scores()) {
					csvRow += buildFscoreColumns(prot.getEvaluationScoreCalculator().getHighestPossibleSimpleGoAnnotationScore());
				}
				if (getSettings().doCalculateAncestryGoF1Scores()) {
					csvRow += buildFscoreColumns(prot.getEvaluationScoreCalculator().getHighestPossibleAncestryGoAnnotationScore());
				}
				if (getSettings().doCalculateSemSimGoF1Scores()) {
					csvRow += buildFscoreColumns(prot.getEvaluationScoreCalculator().getHighestPossibleSemSimGoAnnotationScore());
				}
			}

			// Write row to CSV:
			csvRow += "\n";
			bw.write(csvRow);

			// If AHRD is requested to write out the AHRD-Score of each
			// BlastHit's Description, do so into another file:
			if (getSettings().doWriteHRDScoresToOutput())
				writeHrdScoresOutput(prot);
		}

		// CLEAN UP:
		bw.close();
		if (getSettings().doWriteHRDScoresToOutput())
			this.hrdScoresWriter.close();
	}
	
	public String buildBestBlastHitsHeader() {
		String hdr = "";
		for (String blastDb : getSettings().getBlastDatabases()) {
			if (blastDb != null && !blastDb.equals("")) {
				hdr += ("\tBest-BlastHit-against-" + blastDb + "-Accession-and-Description");
				if (getSettings().isInTrainingMode()) {
					hdr += "\tBest-BlastHit-against-" + blastDb + "-Description-Length\tBest-BlastHit-against-" + blastDb + "-Description-Evaluation-Score";
					if (getSettings().doWriteFscoreDetailsToOutput()) {
						hdr += "\tBest-BlastHit-against-" + blastDb + "-Description-Evaluation-Score-Precision\tBest-BlastHit-against-" + blastDb + "-Description-Evaluation-Score-Recall\tBest-BlastHit-against-" + blastDb + "-Description-Evaluation-Score-Remaining-Uncertainty\tBest-BlastHit-against-" + blastDb + "-Description-Evaluation-Score-Missing-Information";
					}
					if (getSettings().hasGeneOntologyAnnotations() && getSettings().hasGroundTruthGoAnnotations()) {
						hdr += "\tBest-BlastHit-against-" + blastDb + "-GO-Annotations";
						if (getSettings().doCalculateSimpleGoF1Scores()) {
							hdr += "\tBest-BlastHit-against-" + blastDb + "-GO-Annotations-Simple-F-Score";
							if (getSettings().doWriteFscoreDetailsToOutput()) {
								hdr += "\tBest-BlastHit-against-" + blastDb + "-GO-Annotations-Simple-F-Score-Precision";
								hdr += "\tBest-BlastHit-against-" + blastDb + "-GO-Annotations-Simple-F-Score-Recall";
								hdr += "\tBest-BlastHit-against-" + blastDb + "-GO-Annotations-Simple-F-Score-Remaining-Uncertainty";
								hdr += "\tBest-BlastHit-against-" + blastDb + "-GO-Annotations-Simple-F-Score-Missing-Information";
							}
						}
						if (getSettings().doCalculateAncestryGoF1Scores()) {
							hdr += "\tBest-BlastHit-against-" + blastDb + "-GO-Annotations-Ancestry-F-Score";
							if (getSettings().doWriteFscoreDetailsToOutput()) {
								hdr += "\tBest-BlastHit-against-" + blastDb + "-GO-Annotations-Ancestry-F-Score-Precision";
								hdr += "\tBest-BlastHit-against-" + blastDb + "-GO-Annotations-Ancestry-F-Score-Recall";
								hdr += "\tBest-BlastHit-against-" + blastDb + "-GO-Annotations-Ancestry-F-Score-Remaining-Uncertainty";
								hdr += "\tBest-BlastHit-against-" + blastDb + "-GO-Annotations-Ancestry-F-Score-Missing-Information";
							}
						}
						if (getSettings().doCalculateSemSimGoF1Scores()) {
							hdr += "\tBest-BlastHit-against-" + blastDb + "-GO-Annotations-SemSim-F-Score";
							if (getSettings().doWriteFscoreDetailsToOutput()) {
								hdr += "\tBest-BlastHit-against-" + blastDb + "-GO-Annotations-SemSim-F-Score-Precision";
								hdr += "\tBest-BlastHit-against-" + blastDb + "-GO-Annotations-SemSim-F-Score-Recall";
								hdr += "\tBest-BlastHit-against-" + blastDb + "-GO-Annotations-SemSim-F-Score-Remaining-Uncertainty";
								hdr += "\tBest-BlastHit-against-" + blastDb + "-GO-Annotations-SemSim-F-Score-Missing-Information";
							}
						}
					}
				}
			}
		}
		return hdr;
	}

	/**
	 * @param Protein
	 *            prot
	 * @return String - Part of the CSV-Row with the columns Evaluator-Score and
	 *         Ground-Truth-Description.
	 */
	public String buildDescriptionEvaluationColumns(Protein prot) {
		// HEADER:
		// \tHRD-Length\tGround-Truth-Description\tRef-Lenght\tEvaluation-Score\tDiff-to-bestCompetitor
		String csvCells = "";
		// HRD-Length ground truth and AHRD's performance:
		if (prot.getEvaluationScoreCalculator().getEvalutionScore() != null) {
			// HRD-Length ref f1score diff-to-best-competitor:
			csvCells += "\t";
			if (prot.getDescriptionScoreCalculator().getHighestScoringBlastResult() != null) {
				csvCells += prot.getDescriptionScoreCalculator().getHighestScoringBlastResult().getEvaluationTokens().size();
			}
			else {
				csvCells += "0";
			}
			csvCells += "\t" + prot.getEvaluationScoreCalculator().getGroundTruthDescription().getDescription() + "\t"
					+ prot.getEvaluationScoreCalculator().getGroundTruthDescription().getTokens().size() + "\t"
					+ FRMT.format(prot.getEvaluationScoreCalculator().getEvalutionScore().getScore());
			if (getSettings().doWriteFscoreDetailsToOutput()) {
			csvCells += "\t" + FRMT.format(prot.getEvaluationScoreCalculator().getEvalScoreMinBestCompScore()) + "\t"
					+ FRMT.format(prot.getEvaluationScoreCalculator().getEvalutionScore().getPrecision()) + "\t"
					+ FRMT.format(prot.getEvaluationScoreCalculator().getEvalutionScore().getRecall()) + "\t"
					+ FRMT.format(prot.getEvaluationScoreCalculator().getEvalutionScore().getRemainingUncertainty()) + "\t"
					+ FRMT.format(prot.getEvaluationScoreCalculator().getEvalutionScore().getMissingInformation());
			}
		} else {
			csvCells = "\t\t\t\t";
			if (getSettings().doWriteFscoreDetailsToOutput()) {
				csvCells += "\t\t\t\t\t";
			}
		}
		return csvCells;
	}
	
	public String buildTokenSetCell(Protein prot) {
		String tokenSetCell = "\t";
		for (String token : prot.getTokenScoreCalculator().getTokenScores().keySet()) {
			tokenSetCell += "[" + token + "->" + FRMT.format(prot.getTokenScoreCalculator().getTokenScores().get(token))
					+ "]";
		}
		return tokenSetCell;
	}

	/**
	 * Append the following columns to the current CSV-Row:
	 * sum_of_all_token_scores, token_high_score, correction_factor, go_score,
	 * lexical_score, rel_bit_score, desc_line_frequency,
	 * max_desc_line_frequency
	 * 
	 * @note: Printing out these values may slow down the running-time, as these
	 *        values are not necessarily stored in memory.
	 * 
	 * @param prot
	 * @return String - Part of the CSV-Row holding the values of the above
	 *         columns.
	 */
	public String buildDescScoreCells(Protein prot) {
		String csvCells = "";
		// Found a high scoring description?
		if (prot.getDescriptionScoreCalculator().getHighestScoringBlastResult() == null) {
			csvCells = "\t\t\t\t\t";
		} else {
			BlastResult hsbr = prot.getDescriptionScoreCalculator().getHighestScoringBlastResult();
			csvCells += "\t" + FRMT.format(prot.getTokenScoreCalculator().sumOfAllTokenScores(hsbr));
			csvCells += "\t" + FRMT.format(prot.getTokenScoreCalculator().getTokenHighScore());
			csvCells += "\t" + FRMT.format(prot.getLexicalScoreCalculator().correctionFactor(hsbr));
			csvCells += "\t" + FRMT.format(prot.getLexicalScoreCalculator().lexicalScore(hsbr));
			csvCells += "\t" + FRMT.format(prot.getDescriptionScoreCalculator().relativeBlastScore(hsbr));
		}
		return csvCells;
	}

	private String buildFscoreColumns(Fscore f) {
		String csvCells = "\t" + FRMT.format(f.getScore());
		if (getSettings().doWriteFscoreDetailsToOutput()) {
			csvCells += "\t" + FRMT.format(f.getPrecision());
			csvCells += "\t" + FRMT.format(f.getRecall());
			csvCells += "\t" + FRMT.format(f.getRemainingUncertainty());
			csvCells += "\t" + FRMT.format(f.getMissingInformation());
		}
		return csvCells;
	}
	
	private String buildEmptyFscoreColumns() {
		String csvCells =  "\t0";
		if (getSettings().doWriteFscoreDetailsToOutput()) {
			csvCells += "\t0\t0\t1\t1";
		}
		return csvCells;
	}
	
	public String buildBestBlastHitsColumns(Protein prot) {
		String csvRow = "";
		for (String blastDb : getSettings().getBlastDatabases()) {
			if (prot.getEvaluationScoreCalculator().getBestUnchangedBlastResults().get(blastDb) != null) {
				BlastResult bestBr = prot.getEvaluationScoreCalculator().getBestUnchangedBlastResults().get(blastDb);
				csvRow += "\t\"" + bestBr.getAccession() + " " + bestBr.getDescription() + "\"";
				if (getSettings().isInTrainingMode()) {
					csvRow += "\t" + bestBr.getEvaluationTokens().size();
					csvRow += buildFscoreColumns(bestBr.getEvaluationScore());
					if (getSettings().hasGeneOntologyAnnotations() && getSettings().hasGroundTruthGoAnnotations()) {
						csvRow += "\t" + combineGoTermsToString(bestBr.getGoAnnotations());
						if (getSettings().doCalculateSimpleGoF1Scores()) {
							csvRow += buildFscoreColumns(bestBr.getSimpleGoAnnotationScore());
						}
						if (getSettings().doCalculateAncestryGoF1Scores()) {
							csvRow += buildFscoreColumns(bestBr.getAncestryGoAnnotationScore());
						}
						if (getSettings().doCalculateSemSimGoF1Scores()) {
							csvRow += buildFscoreColumns(bestBr.getSemSimGoAnnotationScore());
						}
					}
				}
			} else {
				csvRow += "\t";
				if (getSettings().isInTrainingMode()) {
					csvRow += "\t0";
					csvRow += buildEmptyFscoreColumns();
					if (getSettings().hasGeneOntologyAnnotations() && getSettings().hasGroundTruthGoAnnotations()) {
						csvRow += "\t";
						if (getSettings().doCalculateSimpleGoF1Scores()) {
							csvRow += buildEmptyFscoreColumns();
						}
						if (getSettings().doCalculateAncestryGoF1Scores()) {
							csvRow += buildEmptyFscoreColumns();
						}
						if (getSettings().doCalculateSemSimGoF1Scores()) {
							csvRow += buildEmptyFscoreColumns();
						}
					}
				}
			}
		}
		return csvRow;
	}

	public String buildCompetitorColumns(String competitor, Protein prot) {
		String csvCols = "";
		Map<String, CompetitorAnnotation> compAnnots = prot.getEvaluationScoreCalculator().getCompetitorAnnotations();
		if (compAnnots != null) {
			CompetitorAnnotation annot = compAnnots.get(competitor);
			if (annot != null) {
				csvCols += "\t" + annot.getDescription() + "\t" + annot.getEvaluationTokens().size();
				csvCols += buildFscoreColumns(annot.getEvaluationScore());
				if (getSettings().hasGeneOntologyAnnotations() && getSettings().hasGroundTruthGoAnnotations()) {
					csvCols += "\t" + combineGoTermsToString(annot.getGoAnnotations());
					if (getSettings().doCalculateSimpleGoF1Scores()) {
						csvCols += buildFscoreColumns(annot.getSimpleGoAnnotationScore());
					}
					if (getSettings().doCalculateAncestryGoF1Scores()) {
						csvCols += buildFscoreColumns(annot.getAncestryGoAnnotationScore());
					}
					if (getSettings().doCalculateSemSimGoF1Scores()) {
						csvCols += buildFscoreColumns(annot.getSemSimGoAnnotationScore());
					}
				}
			} else {
				csvCols += buildCompetitorsMissingAnnotationColumns();
			}
		} else {
			csvCols += buildCompetitorsMissingAnnotationColumns();
		}
		return csvCols;
	}

	public String buildCompetitorsMissingAnnotationColumns() {
		String csvCols = "";
		csvCols += "\t\t0";
		csvCols += buildEmptyFscoreColumns();
		if (getSettings().hasGeneOntologyAnnotations() && getSettings().hasGroundTruthGoAnnotations()) {
			csvCols += "\t";
			if (getSettings().doCalculateSimpleGoF1Scores()) {
				csvCols += buildEmptyFscoreColumns();
			}
			if (getSettings().doCalculateAncestryGoF1Scores()) {
				csvCols += buildEmptyFscoreColumns();
			}
			if (getSettings().doCalculateSemSimGoF1Scores()) {
				csvCols += buildEmptyFscoreColumns();
			}
		}
		return csvCols;
	}
	
	public String buildBlastResultWithHighestPossibleDescriptionScoreColumns(Protein prot) {
		String csvCols = "";
		BlastResult br = prot.getEvaluationScoreCalculator().getBlastResultWithHighestPossibleDescriptionScore();
		if (br != null) {
			csvCols += "\t\"" + br.getAccession() + " " + br.getDescription() + "\"";
		} else {
			csvCols += "\t";
		}
		csvCols += buildFscoreColumns(prot.getEvaluationScoreCalculator().getHighestPossibleDescriptionScore());
		return csvCols;
	}

	private String buildGroundTruthGoAnnotationColumns(Protein prot) {
		String goColumns = "\t"
				+ combineGoTermsToString(prot.getEvaluationScoreCalculator().getGroundTruthGoAnnoatations());
		if (getSettings().doCalculateSimpleGoF1Scores()) {
			goColumns += buildFscoreColumns(prot.getEvaluationScoreCalculator().getSimpleGoAnnotationScore());
		}
		if (getSettings().doCalculateAncestryGoF1Scores()) {
			goColumns += buildFscoreColumns(prot.getEvaluationScoreCalculator().getAncestryGoAnnotationScore());
		}
		if (getSettings().doCalculateSemSimGoF1Scores()) {
			goColumns += buildFscoreColumns(prot.getEvaluationScoreCalculator().getSemSimGoAnnotationScore());
		}
		return goColumns;
	}
	
}
