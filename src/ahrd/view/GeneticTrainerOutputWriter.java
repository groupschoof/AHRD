package ahrd.view;

import java.io.IOException;

import ahrd.controller.Parameters;

public class GeneticTrainerOutputWriter extends TrainerOutputWriter {

	public GeneticTrainerOutputWriter() throws IOException {
		super();
	}

	@Override
	public String generateHeader(boolean isFinalOutput) {
		String hdr = "Generation\t";
		if (isFinalOutput)
			hdr += "Average Maximum-Evaluation-Score\t";
		hdr += "Average Evaluation-Score(F-Score)";
		if (!isFinalOutput)
			hdr += "\tDiff-to-last-Generation\tOrigin";
		hdr += "\tAverage True-Positive-Rate\tAverage False-Positive-Rate\tToken-Score-Bit-Score-Weight\tToken-Score-Database-Score-Weight\tToken-Score-Overlap-Score-Weight";
		for (String blastDb : this.sortedBlastDatabases) {
			hdr += "\t" + blastDb + "-Weight";
			hdr += "\t" + blastDb + "-Description-Score-Bit-Score-Weight";
		}
		hdr += "\n";
		return hdr;
	}

	public void writeIterationOutput(int generation, Parameters currentGenerationsBestParameters,
			double diffAvgEvalScoreToLastGeneration, String Origin)
			throws IOException {
		this.pathBufWrtr.write(settingsRow(generation, currentGenerationsBestParameters,
				diffAvgEvalScoreToLastGeneration, Origin));
	}

	public String settingsRow(int generation, Parameters p,
			double diffAvgEvalScoreToLastGeneration, String Origin) {
		String col = generation + "\t"
				+ p.getAvgEvaluationScore() + "\t"
				+ diffAvgEvalScoreToLastGeneration + "\t" + Origin + "\t"
				+ FRMT.format(p.getAvgTruePositivesRate()) + "\t"
				+ FRMT.format(p.getAvgFalsePositivesRate()) + "\t"
				+ FRMT.format(p.getTokenScoreBitScoreWeight()) + "\t"
				+ FRMT.format(p.getTokenScoreDatabaseScoreWeight()) + "\t"
				+ FRMT.format(p.getTokenScoreOverlapScoreWeight());
		for (String blastDb : this.sortedBlastDatabases) {
			col += "\t" + FRMT.format(p.getBlastDbWeight(blastDb));
			col += "\t"
					+ FRMT.format(p.getDescriptionScoreBitScoreWeight(blastDb));
		}
		col += "\n";
		return col;
	}

}
