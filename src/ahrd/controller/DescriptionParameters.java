package ahrd.controller;

import static ahrd.controller.Settings.getSettings;
import static ahrd.controller.Utils.randomMultipleOfOne;
import static ahrd.controller.Utils.randomMultipleOfTen;
import static ahrd.controller.Utils.randomSaveSubtract;
import static ahrd.controller.Utils.roundToNDecimalPlaces;

import java.text.DecimalFormat;
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
 * @author Florian Boecker, Kathrin Klee, Asis Hallab
 */
public class DescriptionParameters extends Parameters implements Cloneable {
	
	@Override
	public int numberOfNonDbParameters() {
		return 3;
	}
	
	/**
	 * 
	 * @param sortedDistinctBlastDatabaseNames
	 * @return
	 */
	@Override
	public Parameters randomParameters(List<String> sortedDistinctBlastDatabaseNames) {
		DescriptionParameters out = new DescriptionParameters();
		Random rand = Utils.random;
		// draw random token score weights
		out.setTokenScoreBitScoreWeight(roundToNDecimalPlaces(rand.nextDouble(), 4));
		out.setTokenScoreDatabaseScoreWeight(roundToNDecimalPlaces(rand.nextDouble(), 4));
		out.setTokenScoreOverlapScoreWeight(roundToNDecimalPlaces(rand.nextDouble(), 4));
		// normalize the randomly chosen weights:
		out.normalizeTokenScoreWeights();
		// Init BlastDbs' Parameters:
		for (String blastDbName : sortedDistinctBlastDatabaseNames) {
			out.setAnnotationScoreBitScoreWeight(blastDbName, randomMultipleOfOne().toString());
			out.setBlastDbWeight(blastDbName, randomMultipleOfTen().toString());
		}
		// Set origin for genetic training output
		out.setOrigin("random");
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
	 * <li>Informative-Token-Threshold</li>
	 * <li>Blast-Database-Weight(different for each Blast-Database)</li>
	 * <li>Description-Score-Bit-Score-Weight (different for each Blast-Database)</li>
	 * </ul>
	 * 
	 * @NOTE: The three <em>Token-Score-Weights</em> <strong>must</strong> sum
	 *        up to 1.
	 * 
	 * @param Double
	 *            diffEvalScoreToLastEvaluatedParams - To increase probability
	 *            to perform 'hill climbing' during optimization, a good
	 *            increase in evaluation score increases likelihood to mutate
	 *            the last mutated parameter again.
	 * @return clone of this instance with one of the above mentioned parameters
	 *         <em>slightly</em> changed.
	 */
	@Override
	public DescriptionParameters neighbour(Double diffEvalScoreToLastEvaluatedParams) {
		DescriptionParameters ngb = this.clone();
		// Randomly decide to mutate the same parameter again, if last mutation
		// resulted in an increase of score:
		Integer randParamToMutate = getLastMutatedParameter();
		if (!(diffEvalScoreToLastEvaluatedParams != null
				&& diffEvalScoreToLastEvaluatedParams > 0.0
				&& randParamToMutate != null && Utils.random.nextDouble() <= pMutateSameParameter(diffEvalScoreToLastEvaluatedParams))) {
			// Do not mutate the same parameter again, but randomly choose one
			// to change:
			randParamToMutate = parameterToMutateRandomIndex();
		}
		// Once a parameter is chosen by its index, mutate it:
		if (randParamToMutate < numberOfNonDbParameters()) {
			// Mutate one of the six parameters independent of the number of Blast-Databases:
			switch(randParamToMutate) {
				case 0: ngb.mutateTokenScoreBitScoreWeight(); break;
				case 1: ngb.mutateTokenScoreDatabaseScoreWeight(); break;
				case 2: ngb.mutateTokenScoreOverlapScoreWeight(); break;
			}
		} else {
			// Mutate a Parameter associated with a Blast-Database:
			int indOfBlastDbToMutate = randParamToMutate - numberOfNonDbParameters();
			int blastDbIndex = (new Double(
					Math.floor(indOfBlastDbToMutate / 2.0))).intValue();
			String blastDbToMutate = getSettings().getSortedBlastDatabases()
					.get(blastDbIndex);
			boolean mutateWeight = (indOfBlastDbToMutate % 2 == 0);
			if (mutateWeight)
				ngb.mutateBlastDatabaseWeight(blastDbToMutate);
			else
				ngb.mutateAnnotationScoreBitScoreWeight(blastDbToMutate);
		}
		// Remember what made the neighbor different from its parent:
		ngb.setLastMutatedParameter(randParamToMutate);
		// Reset average evaluation score
		ngb.setAvgEvaluationScore(null);
		ngb.setAvgPrecision(null);
		ngb.setAvgRecall(null);
		// Set origin for genetic training output
		ngb.setOrigin("mutation");
		return ngb;
	}

	/**
	 * Creates an offspring with a random recombination of the current parameters and a given parameter set.
	 * 
	 * @NOTE: The three <em>Token-Score-Weights</em> are normalized to sum up to 1.
	 * 
	 * @param partner - The Parameters to recombine the current ones with 
	 * 
	 * @return The random offspring of the current and given Parameters  
	 */
	@Override
	public Parameters recombine(Parameters partner) {
		DescriptionParameters offspring = this.clone();
		Random rand = Utils.random;
		if(rand.nextBoolean())
			offspring.setTokenScoreBitScoreWeight(partner.getTokenScoreBitScoreWeight());
		if(rand.nextBoolean())
			offspring.setTokenScoreDatabaseScoreWeight(partner.getTokenScoreDatabaseScoreWeight());
		if(rand.nextBoolean())
			offspring.setTokenScoreOverlapScoreWeight(partner.getTokenScoreOverlapScoreWeight());
		for (String blastDbName : getSettings().getSortedBlastDatabases()) {
			if(rand.nextBoolean())
				offspring.setAnnotationScoreBitScoreWeight(blastDbName, partner.getAnnotationScoreBitScoreWeight(blastDbName).toString());
			if(rand.nextBoolean())
				offspring.setBlastDbWeight(blastDbName, partner.getBlastDbWeight(blastDbName).toString());
		}
		offspring.normalizeTokenScoreWeights();
		offspring.setAvgEvaluationScore(null);
		offspring.setAvgPrecision(null);
		offspring.setAvgRecall(null);
		// Set origin for genetic training
		offspring.setOrigin("recombination");
		return offspring;
	}

	/**
	 * Returns a clone of this instance.
	 */
	@Override
	public DescriptionParameters clone() {
		DescriptionParameters clone = (DescriptionParameters) super.clone();
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
		if (!(eql instanceof DescriptionParameters))
			return false;
		// We are dealing with an Instance of Parameters:
		boolean areBlastParamsEqual = true;
		for (String blastDb : getBlastDbParameters().keySet()) {
			for (String iterKey : getParametersOfBlastDb(blastDb).keySet()) {
				areBlastParamsEqual = areBlastParamsEqual
						&& getParametersOfBlastDb(blastDb).get(iterKey).equals(
								((DescriptionParameters) eql).getParametersOfBlastDb(
										blastDb).get(iterKey));
			}
		}
		return areBlastParamsEqual
				&& ((DescriptionParameters) eql).getTokenScoreBitScoreWeight().equals(
						this.getTokenScoreBitScoreWeight())
				&& ((DescriptionParameters) eql).getTokenScoreDatabaseScoreWeight()
						.equals(this.getTokenScoreDatabaseScoreWeight())
				&& ((DescriptionParameters) eql).getTokenScoreOverlapScoreWeight().equals(this.getTokenScoreOverlapScoreWeight());
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
				+ getTokenScoreOverlapScoreWeight();
		return hashSrc.hashCode();
	}
	
	@Override
	public String buildHeaderForOutput(String separator) {
		String hdr = "Token-Score-Bit-Score-Weight" + separator
				+ "Token-Score-Database-Score-Weight" + separator
				+ "Token-Score-Overlap-Score-Weight";
		for (String blastDb : blastDbParameters.keySet()) {
			hdr += separator + blastDb + "-Weight";
			hdr += separator + blastDb + "-Description-Score-Bit-Score-Weight";
		}
		return hdr;
	}
	
	@Override
	public String formatForOutput(DecimalFormat frmt, String separator) {
		String out = frmt.format(getTokenScoreBitScoreWeight()) + separator
				+ frmt.format(getTokenScoreDatabaseScoreWeight()) + separator
				+ frmt.format(getTokenScoreOverlapScoreWeight());
		for (String blastDb : blastDbParameters.keySet()) {
			out += separator + frmt.format(getBlastDbWeight(blastDb));
			out += separator + frmt.format(getAnnotationScoreBitScoreWeight(blastDb));
		}
		return out;
	}

}