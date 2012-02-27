package ahrd.controller;

import static ahrd.controller.Utils.roundToNDecimalPlaces;
import static ahrd.controller.Utils.randomMultipleOfOneTenth;
import static ahrd.controller.Utils.randomMultipleOfTen;
import static ahrd.controller.Utils.randomSaveSubtract;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

/**
 * The following parameters are those subject to optimization. They are stored
 * wrapped in a distinct class from Settings in order to enable random
 * generation and scoring of these parameters.
 * 
 * @author Kathrin Klee, Asis Hallab
 */
public class Parameters implements Cloneable {

	private Double tokenScoreBitScoreWeight;
	private Double tokenScoreDatabaseScoreWeight;
	private Double tokenScoreOverlapScoreWeight;
	private Double descriptionScorePatternFactorWeight;
	private Map<String, Map<String, String>> blastDbParameters = new HashMap<String, Map<String, String>>();
	/**
	 * If we test different settings in the parameter-space, remember the
	 * average evaluation-score (objective-function).
	 */
	private Double avgEvaluationScore;
	/**
	 * If we test different settings in the parameter-space, remember the
	 * average True-Positives-Rate (TPR).
	 */
	private Double avgTruePositivesRate;
	/**
	 * If we test different settings in the parameter-space, remember the
	 * average False-Positives-Rate (FPR).
	 */
	private Double avgFalsePositivesRate;

	public static Parameters randomParameters(
			List<String> sortedDistinctBlastDatabaseNames) {
		Parameters out = new Parameters();
		out.setTokenScoreBitScoreWeight(randomMultipleOfOneTenth());
		out.setTokenScoreDatabaseScoreWeight(randomMultipleOfOneTenth());
		out.setTokenScoreOverlapScoreWeight(randomMultipleOfOneTenth());
		// normalize the randomly chosen weights:
		out.normalizeTokenScoreWeights();
		out.setDescriptionScorePatternFactorWeight(randomMultipleOfOneTenth());
		// Init BlastDbs' Parameters:
		for (String blastDbName : sortedDistinctBlastDatabaseNames) {
			out.setDescriptionScoreBitScoreWeight(blastDbName,
					randomMultipleOfOneTenth().toString());
			out.setBlastDbWeight(blastDbName, randomMultipleOfTen().toString());
		}
		return out;
	}

	/**
	 * Clones this instance and changes one of the following <em>six</em>
	 * fields, AHRD-parameters, in order to calculate AHRD's performance with
	 * different parameters:
	 * <ul>
	 * <li>Token-Score-Bit-Score-Weight</li>
	 * <li>Token-Score-Database-Score-Weight</li>
	 * <li>Token-Score-Overlap-Score-Weight</li>
	 * <li>Blast-Database-Weight</li>
	 * <li>Description-Score-Bit-Score-Weight (different for each
	 * Blast-Database)</li>
	 * <li>Description-Score-Relative-Description-Frequency-Weight</li>
	 * </ul>
	 * 
	 * @NOTE: The three <em>Token-Score-Weights</em> <strong>must</strong> sum
	 *        up to 1.
	 * 
	 * @return clone of this instance with one of the above mentioned parameters
	 *         <em>slightly</em> changed.
	 */
	public Parameters neighbour() {
		Parameters ngb = this.clone();
		// Randomly choose a parameter to change:
		Random rand = Utils.random;
		int randParamInd = rand.nextInt(6);
		switch (randParamInd) {
		case 0:
			ngb.mutateTokenScoreBitScoreWeight();
			break;
		case 1:
			ngb.mutateTokenScoreDatabaseScoreWeight();
			break;
		case 2:
			ngb.mutateTokenScoreOverlapScoreWeight();
			break;
		case 3:
			ngb.mutateBlastDatabaseWeight();
			break;
		case 4:
			ngb.mutateDescriptionScoreBitScoreWeight();
			break;
		case 5:
			ngb.mutateDescriptionScorePatternFactorWeight();
			break;
		}
		return ngb;
	}

	public String randomBlastDatabaseName() {
		Random rand = Utils.random;
		int randBlastDbInd = rand.nextInt(getBlastDatabases().size());
		List<String> blastDbNamesList = new ArrayList<String>(
				getBlastDatabases());
		return blastDbNamesList.get(randBlastDbInd);
	}

	public void mutateBlastDatabaseWeight() {
		String blastDb = randomBlastDatabaseName();
		Integer bdbw = getBlastDbWeight(blastDb);
		if (randomSaveSubtract(bdbw, Settings.BLAST_DB_WEIGHT_MUTATOR_SEED))
			bdbw -= Settings.BLAST_DB_WEIGHT_MUTATOR_SEED;
		else
			bdbw += Settings.BLAST_DB_WEIGHT_MUTATOR_SEED;

		setBlastDbWeight(blastDb, bdbw.toString());
	}

	public void mutateDescriptionScoreBitScoreWeight() {
		String blastDb = randomBlastDatabaseName();
		Double bsw = getDescriptionScoreBitScoreWeight(blastDb);
		if (randomSaveSubtract(bsw, Settings.PERCENTAGE_MUTATOR_SEED))
			bsw -= Settings.PERCENTAGE_MUTATOR_SEED;
		else
			bsw += Settings.PERCENTAGE_MUTATOR_SEED;
		setDescriptionScoreBitScoreWeight(blastDb, bsw.toString());
	}

	public void mutateDescriptionScorePatternFactorWeight() {
		double pfw = getDescriptionScorePatternFactorWeight();
		if (randomSaveSubtract(pfw, Settings.PERCENTAGE_MUTATOR_SEED))
			pfw -= Settings.PERCENTAGE_MUTATOR_SEED;
		else
			pfw += Settings.PERCENTAGE_MUTATOR_SEED;
		setDescriptionScorePatternFactorWeight(pfw);
	}

	/**
	 * Normalizes the three weights appearing in the Token-Score-Formula, so
	 * they sum up to 1.0
	 */
	public void normalizeTokenScoreWeights() {
		double s = roundToNDecimalPlaces(getTokenScoreBitScoreWeight()
				+ getTokenScoreDatabaseScoreWeight()
				+ getTokenScoreOverlapScoreWeight(), 4);
		setTokenScoreBitScoreWeight(roundToNDecimalPlaces(
				getTokenScoreBitScoreWeight() / s, 4));
		setTokenScoreDatabaseScoreWeight(roundToNDecimalPlaces(
				getTokenScoreDatabaseScoreWeight() / s, 4));
		setTokenScoreOverlapScoreWeight(roundToNDecimalPlaces(
				getTokenScoreOverlapScoreWeight() / s, 4));
	}

	/**
	 * Diminishes or increases Token-Score-Bit-Score-Weight by
	 * PERCENTAGE_MUTATOR_SEED and normalizes the other two weights in the
	 * Token-Score-Formula.
	 */
	public void mutateTokenScoreBitScoreWeight() {
		Double bsw = getTokenScoreBitScoreWeight();
		if (randomSaveSubtract(bsw, Settings.PERCENTAGE_MUTATOR_SEED))
			bsw = bsw - Settings.PERCENTAGE_MUTATOR_SEED;
		else
			bsw = bsw + Settings.PERCENTAGE_MUTATOR_SEED;
		setTokenScoreBitScoreWeight(bsw);
		// normalize:
		normalizeTokenScoreWeights();
	}

	/**
	 * Diminishes or increases Token-Score-Database-Score-Weight by
	 * PERCENTAGE_MUTATOR_SEED and normalizes the other two weights in the
	 * Token-Score-Formula.
	 */
	public void mutateTokenScoreDatabaseScoreWeight() {
		Double dbsw = getTokenScoreDatabaseScoreWeight();
		if (randomSaveSubtract(dbsw, Settings.PERCENTAGE_MUTATOR_SEED))
			dbsw = dbsw - Settings.PERCENTAGE_MUTATOR_SEED;
		else
			dbsw = dbsw + Settings.PERCENTAGE_MUTATOR_SEED;
		setTokenScoreDatabaseScoreWeight(dbsw);
		// normalize:
		normalizeTokenScoreWeights();
	}

	/**
	 * Diminishes or increases Token-Score-Overlap-Score-Weight by
	 * PERCENTAGE_MUTATOR_SEED and normalizes the other two weights in the
	 * Token-Score-Formula.
	 */
	public void mutateTokenScoreOverlapScoreWeight() {
		Double osw = getTokenScoreOverlapScoreWeight();
		if (randomSaveSubtract(osw, Settings.PERCENTAGE_MUTATOR_SEED))
			osw = osw - Settings.PERCENTAGE_MUTATOR_SEED;
		else
			osw = osw + Settings.PERCENTAGE_MUTATOR_SEED;
		setTokenScoreOverlapScoreWeight(osw);
		// normalize:
		normalizeTokenScoreWeights();
	}

	/**
	 * Returns a clone of this instance.
	 */
	public Parameters clone() {
		Parameters clone;
		try {
			clone = (Parameters) super.clone();
		} catch (CloneNotSupportedException e) {
			e.printStackTrace(System.err);
			return null;
		}
		// Clone the Blast-Database-Parameters-Map and Values:
		Map<String, Map<String, String>> blastDbSettings = new HashMap<String, Map<String, String>>();
		for (String blastDb : getBlastDbParameters().keySet()) {
			blastDbSettings.put(blastDb, new HashMap<String, String>());
			for (String iterKey : getParametersOfBlastDb(blastDb).keySet()) {
				blastDbSettings.get(blastDb)
						.put(new String(iterKey),
								new String(getParametersOfBlastDb(blastDb).get(
										iterKey)));
			}
		}
		clone.blastDbParameters = blastDbSettings;
		return clone;
	}

	@Override
	public boolean equals(Object eql) {
		if (!(eql instanceof Parameters))
			return false;
		// We are dealing with an Instance of Parameters:
		boolean areBlastParamsEqual = true;
		for (String blastDb : getBlastDbParameters().keySet()) {
			for (String iterKey : getParametersOfBlastDb(blastDb).keySet()) {
				areBlastParamsEqual = areBlastParamsEqual
						&& getParametersOfBlastDb(blastDb).get(iterKey).equals(
								((Parameters) eql).getParametersOfBlastDb(
										blastDb).get(iterKey));
			}
		}
		return areBlastParamsEqual
				&& ((Parameters) eql).getTokenScoreBitScoreWeight().equals(
						this.getTokenScoreBitScoreWeight())
				&& ((Parameters) eql).getTokenScoreDatabaseScoreWeight()
						.equals(this.getTokenScoreDatabaseScoreWeight())
				&& ((Parameters) eql).getTokenScoreOverlapScoreWeight().equals(
						this.getTokenScoreOverlapScoreWeight())
				&& ((Parameters) eql).getDescriptionScorePatternFactorWeight()
						.equals(this.getDescriptionScorePatternFactorWeight());
	}

	@Override
	public int hashCode() {
		String hashSrc = "";
		for (String blastDb : getBlastDbParameters().keySet()) {
			for (String iterKey : getParametersOfBlastDb(blastDb).keySet()) {
				hashSrc += getParametersOfBlastDb(blastDb).get(iterKey);
			}
		}
		hashSrc += getTokenScoreBitScoreWeight()
				+ getTokenScoreDatabaseScoreWeight()
				+ getTokenScoreOverlapScoreWeight()
				+ getDescriptionScorePatternFactorWeight();
		return hashSrc.hashCode();
	}

	/**
	 * @return Set<String> the names of the blast-databases used in the current
	 *         AHRD-Run.
	 */
	public Set<String> getBlastDatabases() {
		return getBlastDbParameters().keySet();
	}

	protected Map<String, String> getParametersOfBlastDb(String blastDbName) {
		Map<String, String> out = getBlastDbParameters().get(blastDbName);
		// Init new on first request:
		if (out == null) {
			out = new HashMap<String, String>();
			getBlastDbParameters().put(blastDbName, out);
		}
		return out;
	}

	public Integer getBlastDbWeight(String blastDatabaseName) {
		return Integer.parseInt(getParametersOfBlastDb(blastDatabaseName).get(
				Settings.BLAST_DB_WEIGHT_KEY));
	}

	public void setBlastDbWeight(String blastDatabaseName, String bdbw) {
		getParametersOfBlastDb(blastDatabaseName).put(
				Settings.BLAST_DB_WEIGHT_KEY, bdbw);
	}

	public Double getDescriptionScoreBitScoreWeight(String blastDatabaseName) {
		return Double.parseDouble(getParametersOfBlastDb(blastDatabaseName)
				.get(Settings.DESCRIPTION_SCORE_BIT_SCORE_WEIGHT));
	}

	public void setDescriptionScoreBitScoreWeight(String blastDatabaseName,
			String dsbsw) {
		getParametersOfBlastDb(blastDatabaseName).put(
				Settings.DESCRIPTION_SCORE_BIT_SCORE_WEIGHT, dsbsw);
	}

	public Double getTokenScoreBitScoreWeight() {
		return tokenScoreBitScoreWeight;
	}

	public void setTokenScoreBitScoreWeight(Double tokenScoreBitScoreWeight) {
		this.tokenScoreBitScoreWeight = tokenScoreBitScoreWeight;
	}

	public Double getTokenScoreDatabaseScoreWeight() {
		return tokenScoreDatabaseScoreWeight;
	}

	public void setTokenScoreDatabaseScoreWeight(
			Double tokenScoreDatabaseScoreWeight) {
		this.tokenScoreDatabaseScoreWeight = tokenScoreDatabaseScoreWeight;
	}

	public Double getTokenScoreOverlapScoreWeight() {
		return tokenScoreOverlapScoreWeight;
	}

	public void setTokenScoreOverlapScoreWeight(
			Double tokenScoreOverlapScoreWeight) {
		this.tokenScoreOverlapScoreWeight = tokenScoreOverlapScoreWeight;
	}

	public Double getDescriptionScorePatternFactorWeight() {
		return descriptionScorePatternFactorWeight;
	}

	public void setDescriptionScorePatternFactorWeight(
			Double descriptionScorePatternFactorWeight) {
		this.descriptionScorePatternFactorWeight = descriptionScorePatternFactorWeight;
	}

	public Double getAvgEvaluationScore() {
		return avgEvaluationScore;
	}

	public void setAvgEvaluationScore(Double avgEvaluationScore) {
		this.avgEvaluationScore = avgEvaluationScore;
	}

	public Map<String, Map<String, String>> getBlastDbParameters() {
		return blastDbParameters;
	}

	public Double getAvgTruePositivesRate() {
		return avgTruePositivesRate;
	}

	public void setAvgTruePositivesRate(Double avgTruePositivesRate) {
		this.avgTruePositivesRate = avgTruePositivesRate;
	}

	public Double getAvgFalsePositivesRate() {
		return avgFalsePositivesRate;
	}

	public void setAvgFalsePositivesRate(Double avgFalsePositivesRate) {
		this.avgFalsePositivesRate = avgFalsePositivesRate;
	}
}
