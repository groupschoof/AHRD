package ahrd.controller;

import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.esotericsoftware.yamlbeans.YamlException;
import com.esotericsoftware.yamlbeans.YamlReader;

public class TrainerBatcher extends Batcher {

	/**
	 * Generate the following number of batches as different start-positions for
	 * the simulated annealing in Parameter-Space.
	 */
	private Integer noOfBatches = 1000;
	private Set<Parameters> distinctStartPositionsInParameterSpace = new HashSet<Parameters>();

	public TrainerBatcher(Map<String, Object> input) {
		super(input);
		// Initialize number of Start-Positions in Parameter-Space to generate
		// Parameter-Sets for:
		if (getInput().get(Settings.NO_START_POSITIONS_IN_PARAM_SPACE) != null)
			this.noOfBatches = Integer.parseInt((String) input
					.get(Settings.NO_START_POSITIONS_IN_PARAM_SPACE));
	}

	@SuppressWarnings("unchecked")
	public static void main(String[] args) throws YamlException, IOException {
		YamlReader reader = new YamlReader(new FileReader(args[0]));
		TrainerBatcher trainerBatcher = new TrainerBatcher(
				(Map<String, Object>) reader.read());
		trainerBatcher.batch();
		trainerBatcher.writeOutput();
		// Log
		System.out
				.println("For each Batch one AHRD input-file has been written into: '"
						+ trainerBatcher.getInput().get(BATCH_YMLS_DIR_KEY)
								.toString() + "'.");
		System.out
				.println("Created shell-script to start AHRD on all Batches in parallel: '"
						+ trainerBatcher.getInput().get(SHELL_SCRIPT_KEY)
						+ "'.");
	}

	public void batch() {
		for (Integer batchName = 0; batchName < this.noOfBatches; batchName++) {
			getOutput().add(
					generateYml(batchName.toString() + ".yml",
							generateDistinctRandomParameters()));
		}
	}

	@SuppressWarnings("unchecked")
	private List<String> sortedDistinctBlastDatabaseNames() {
		List<String> blastDbs = new ArrayList<String>(
				((Map<String, Object>) getInput().get(Settings.BLAST_DBS_KEY))
						.keySet());
		Collections.sort(blastDbs);
		return blastDbs;
	}

	private Parameters generateDistinctRandomParameters() {
		boolean gotPairwiseDistinctRandParams = false;
		Parameters params = null;
		while (!gotPairwiseDistinctRandParams) {
			params = Parameters
					.randomParameters(sortedDistinctBlastDatabaseNames());
			gotPairwiseDistinctRandParams = this.distinctStartPositionsInParameterSpace
					.add(params);
		}
		return params;
	}

	/**
	 * Argument batchName is expected to be the complete file-name!
	 */
	@SuppressWarnings("unchecked")
	public Map<String, Object> generateYml(String batchName,
			Parameters parameters) {
		Map<String, Object> batchYml = new HashMap<String, Object>();

		// Store batchName:
		batchYml.put(BATCH_NAME_KEY, batchName);

		// Set proteins_fasta:
		batchYml.put(Settings.PROTEINS_FASTA_KEY,
				getInput().get(Settings.PROTEINS_FASTA_KEY));

		// Set references:
		batchYml.put(Settings.REFERENCES_FASTA_KEY,
				getInput().get(Settings.REFERENCES_FASTA_KEY));

		// Set blast2go-results, if given:
		if (getInput().containsKey(BLAST_2_GO_RESULTS_DIR_KEY)) {
			batchYml.put(Settings.BLAST_2_GO_ANNOT_FILE_KEY,
					getInput().get(Settings.BLAST_2_GO_ANNOT_FILE_KEY));
		}

		// Store F-Score-Beta-Parameter, if given:
		if (getInput().containsKey(Settings.F_MEASURE_BETA_PARAM_KEY)) {
			batchYml.put(Settings.F_MEASURE_BETA_PARAM_KEY,
					getInput().get(Settings.F_MEASURE_BETA_PARAM_KEY));
		}

		// Put weight-parameters:
		batchYml.put(Settings.TOKEN_SCORE_BIT_SCORE_WEIGHT,
				parameters.getTokenScoreBitScoreWeight());

		batchYml.put(Settings.TOKEN_SCORE_DATABASE_SCORE_WEIGHT,
				parameters.getTokenScoreDatabaseScoreWeight());

		batchYml.put(Settings.TOKEN_SCORE_OVERLAP_SCORE_WEIGHT,
				parameters.getTokenScoreOverlapScoreWeight());

		batchYml.put(
				Settings.DESCRIPTION_SCORE_RELATIVE_DESCIPTION_FREQUENCY_WEIGHT,
				parameters.getDescriptionScorePatternFactorWeight());

		// Put data for each Blast-Db into Yml-Hash
		Map<String, Object> inputBlastDbs = (Map<String, Object>) ((HashMap<String, Object>) getInput()
				.get(Settings.BLAST_DBS_KEY)).clone();
		for (String blastDbName : inputBlastDbs.keySet()) {
			// read input:
			Map<String, String> inputBlastDb = (Map<String, String>) ((HashMap<String, String>) inputBlastDbs
					.get(blastDbName)).clone();

			// Populate from argument parameters:
			inputBlastDb.put(Settings.BLAST_DB_WEIGHT_KEY, parameters
					.getBlastDbWeight(blastDbName).toString());
			inputBlastDb.put(Settings.DESCRIPTION_SCORE_BIT_SCORE_WEIGHT,
					parameters.getDescriptionScoreBitScoreWeight(blastDbName)
							.toString());

			// reference to current Blast-Database-Parameters:
			inputBlastDbs.put(blastDbName, inputBlastDb);
		}
		// Store now extended Blast-Database-Parameters
		batchYml.put(Settings.BLAST_DBS_KEY, inputBlastDbs);

		// Interpro-Data:
		batchYml.put(Settings.INTERPRO_DATABASE_KEY,
				getInput().get(Settings.INTERPRO_DATABASE_KEY));

		batchYml.put(Settings.INTERPRO_RESULT_KEY,
				getInput().get(Settings.INTERPRO_RESULT_KEY));

		// Gene-Ontology-Result:
		batchYml.put(Settings.GENE_ONTOLOGY_RESULT_KEY,
				getInput().get(Settings.GENE_ONTOLOGY_RESULT_KEY));

		// Output-File and Log of simulated annealing's path through parameter-
		// and score-space:
		String outputDir = getInput().get(OUTPUT_DIR_KEY).toString();
		if (!outputDir.endsWith("/"))
			outputDir += "/";
		String batchFileName = batchName.replaceAll("\\.\\S+$", "");
		String outputFile = outputDir + batchFileName + "_ahrd_trainer_out.csv";
		batchYml.put(Settings.OUTPUT_KEY, outputFile);
		String simAnnealPathLog = outputDir + batchFileName
				+ "_sim_anneal_path_log.csv";
		batchYml.put(Settings.SIMULATED_ANNEALING_PATH_LOG_KEY,
				simAnnealPathLog);

		// Append best BlastHits to output?:
		String appendBestBlastHitsToOutput = (String) getInput().get(
				Settings.WRITE_BEST_BLAST_HITS_TO_OUTPUT);

		if (appendBestBlastHitsToOutput != null
				&& new Boolean(appendBestBlastHitsToOutput))
			batchYml.put(Settings.WRITE_BEST_BLAST_HITS_TO_OUTPUT, true);

		// Append the set of Tokens?
		String appendTokenSet = (String) getInput().get(
				Settings.WRITE_TOKEN_SET_TO_OUTPUT);

		if (appendTokenSet != null && new Boolean(appendTokenSet)) {
			batchYml.put(Settings.WRITE_TOKEN_SET_TO_OUTPUT, true);
		}

		// Append Description-related scores?
		String appendDescScores = (String) getInput().get(
				Settings.WRITE_SCORES_TO_OUTPUT);

		if (appendDescScores != null && new Boolean(appendDescScores)) {
			batchYml.put(Settings.WRITE_SCORES_TO_OUTPUT, true);
		}

		// Pass on simulated annealing parameters, if any are given:
		if (getInput().get(Settings.TEMPERATURE_KEY) != null)
			batchYml.put(Settings.TEMPERATURE_KEY,
					(String) getInput().get(Settings.TEMPERATURE_KEY));
		if (getInput().get(Settings.COOL_DOWN_BY_KEY) != null)
			batchYml.put(Settings.COOL_DOWN_BY_KEY,
					(String) getInput().get(Settings.COOL_DOWN_BY_KEY));
		if (getInput()
				.get(Settings.OPTIMIZATION_ACCEPTANCE_PROBABILITY_SCALING_FACTOR_KEY) != null)
			batchYml.put(
					Settings.OPTIMIZATION_ACCEPTANCE_PROBABILITY_SCALING_FACTOR_KEY,
					(String) getInput()
							.get(Settings.OPTIMIZATION_ACCEPTANCE_PROBABILITY_SCALING_FACTOR_KEY));
		if (getInput().get(Settings.MUTATOR_MEAN_KEY) != null)
			batchYml.put(Settings.MUTATOR_MEAN_KEY,
					(String) getInput().get(Settings.MUTATOR_MEAN_KEY));
		if (getInput().get(Settings.MUTATOR_DEVIATION_KEY) != null)
			batchYml.put(Settings.MUTATOR_DEVIATION_KEY, (String) getInput()
					.get(Settings.MUTATOR_DEVIATION_KEY));

		// Path to Batch's YML-File:
		batchYml.put(PATH_TO_BATCH_YML_KEY, generatePathToBatchYml(batchName));

		// Remember already tested Parameter-Sets? Remember path through
		// simulated annealing?
		if (getInput().get(Settings.REMEMBER_SIMULATED_ANNEALING_PATH_KEY) != null
				&& Boolean.parseBoolean(getInput().get(
						Settings.REMEMBER_SIMULATED_ANNEALING_PATH_KEY)
						.toString()))
			batchYml.put(Settings.REMEMBER_SIMULATED_ANNEALING_PATH_KEY, true);

		// Pass on boolean parameter find_highest_possible_evaluation_score, if
		// given and set to true:
		if (getInput().get(Settings.FIND_HIGHEST_POSSIBLE_EVALUATION_SCORE_KEY) != null
				&& Boolean.parseBoolean(getInput().get(
						Settings.FIND_HIGHEST_POSSIBLE_EVALUATION_SCORE_KEY)
						.toString()))
			batchYml.put(Settings.FIND_HIGHEST_POSSIBLE_EVALUATION_SCORE_KEY,
					true);

		return batchYml;
	}
}
