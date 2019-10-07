package ahrd.view;

import java.io.IOException;

import ahrd.controller.Parameters;
import ahrd.controller.Settings;

public class GeneticTrainerOutputWriter extends TrainerOutputWriter {

	public GeneticTrainerOutputWriter() throws IOException {
		super();
	}

	@Override
	public String generateHeader(boolean isFinalOutput) {
		String hdr = "Generation\t";
		if (isFinalOutput)
			hdr += "Average Maximum-Evaluation-Score\t";
		hdr += "Average Training-Score(F-Score)";
		if (!isFinalOutput)
			hdr += "\tDiff-to-last-Generation\tOrigin";
		hdr += "\tToken-Score-Bit-Score-Weight\tToken-Score-Database-Score-Weight\tToken-Score-Overlap-Score-Weight\tGo-Term-Score-Information-Content-Weight\tInformative-Token-Threshold";
		for (String blastDb : this.sortedBlastDatabases) {
			hdr += "\t" + blastDb + "-Weight";
			hdr += "\t" + blastDb + "-Description-Score-Bit-Score-Weight";
		}
		hdr += "\tDescription-Score-Threshold";
		hdr += "\n";
		return hdr;
	}

	public void writeGeneticIterationOutput(int generation, Parameters currentGenerationsBestParameters,
			double diffAvgEvalScoreToLastGeneration, String Origin)
			throws IOException {
		this.pathBufWrtr.write(geneticSettingsRow(generation, currentGenerationsBestParameters,
				diffAvgEvalScoreToLastGeneration, Origin));
		this.pathBufWrtr.flush();
	}

	public String geneticSettingsRow(int generation, Parameters p,
			double diffAvgEvalScoreToLastGeneration, String Origin) {
		String col = generation + "\t"
				+ p.getAvgTrainingScore() + "\t"
				+ diffAvgEvalScoreToLastGeneration + "\t" + Origin + "\t"
				+ FRMT.format(p.getTokenScoreBitScoreWeight()) + "\t"
				+ FRMT.format(p.getTokenScoreDatabaseScoreWeight()) + "\t"
				+ FRMT.format(p.getTokenScoreOverlapScoreWeight()) + "\t"
				+ FRMT.format(p.getGoTermScoreInformationContentWeight()) + "\t"
				+ FRMT.format(p.getInformativeTokenThreshold());
		for (String blastDb : this.sortedBlastDatabases) {
			col += "\t" + FRMT.format(p.getBlastDbWeight(blastDb));
			col += "\t"
					+ FRMT.format(p.getDescriptionScoreBitScoreWeight(blastDb));
		}
		col += "\t" + FRMT.format(p.getDescriptionScoreThreshold());
		col += "\n";
		return col;
	}
	
	@Override
	public String finalSettingsRow(Settings s, Integer sFoundInGeneration,
			Double avgMaxEvalScore) {
		String col = sFoundInGeneration + "\t" + avgMaxEvalScore + "\t"
				+ s.getAvgTrainingScore() + "\t"
				+ s.getTokenScoreBitScoreWeight() + "\t"
				+ s.getTokenScoreDatabaseScoreWeight() + "\t"
				+ s.getTokenScoreOverlapScoreWeight() + "\t"
				+ s.getGoTermScoreInformationContentWeight() + "\t"
				+ s.getInformativeTokenThreshold();
		for (String blastDb : this.sortedBlastDatabases) {
			col += "\t" + s.getBlastDbWeight(blastDb);
			col += "\t" + s.getDescriptionScoreBitScoreWeight(blastDb);
		}
		col += "\t" + s.getDescriptionScoreThreshold();
		col += "\n";
		return col;
	}

}
