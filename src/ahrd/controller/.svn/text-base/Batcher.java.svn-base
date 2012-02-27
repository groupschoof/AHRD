package ahrd.controller;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Pattern;

import com.esotericsoftware.yamlbeans.YamlException;
import com.esotericsoftware.yamlbeans.YamlReader;
import com.esotericsoftware.yamlbeans.YamlWriter;

public class Batcher {

	public static final String AHRD_CALL_KEY = "ahrd_call";
	public static final String AHRD_CALL_BATCH = "#batch#";
	public static final String AHRD_CALL_BATCH_NAME = "#batch_name#";
	public static final String PROTEINS_DIR_KEY = "proteins_dir";
	public static final String REFERENCES_DIR_KEY = "references_dir";
	public static final String BLAST_2_GO_RESULTS_DIR_KEY = "blast2go_dir";
	public static final String BLAST_RESULTS_DIR_KEY = "dir";
	public static final String INTERPRO_RESULTS_DIR_KEY = "interpro_results_dir";
	public static final String INTERPRO_RESULTS_FILE_KEY = "interpro_results_file";
	public static final String GENE_ONTOLOGY_RESULTS_DIR_KEY = "gene_ontology_results_dir";
	public static final String GENE_ONTOLOGY_RESULTS_FILE_KEY = "gene_ontology_results_file";
	public static final String BATCH_YMLS_DIR_KEY = "batch_ymls_dir";
	public static final String OUTPUT_DIR_KEY = "output_dir";
	public static final String SHELL_SCRIPT_KEY = "shell_script";
	public static final String PATH_TO_BATCH_YML_KEY = "path_to_batch_yml";
	public static final String BATCH_NAME_KEY = "batch_name";

	/**
<<<<<<< .mine
	 * Batch-Names e.g. 'batch1' are expected to be just the part of a file-name
	 * before the file-extension, in order to distinguish batch1.fa from
	 * batch11.fa we append "\\." to the regular-expression used to match. We
	 * prepend "\\D" to the regular expression, if and only if the name to match
	 * is a pure number. This is done in order to distinguish 'sprot_1.pairwise'
	 * from 'sprot_11.pairwise', where the Batch-Name to match is just '1'.
=======
	 * Batch-Names are expected to be just the part of a file-name before the
	 * file-extension, in order to distinguish batch1.fa from batch11.fa we
	 * append "\\." to the regular-expression used to match. We prepend "\\D" to
	 * the regular expression, if and only if the name to match is a pure
	 * number. This is done in order to distinguish 'sprot_1.pairwise' from
	 * 'sprot_11.pairwise', where the Batch-Name to match is just '1'.
>>>>>>> .r2783
	 * 
	 * @author hallab, klee
	 */
	public static class BatchFilenameFilter implements FilenameFilter {

		private Pattern matchName;

		public BatchFilenameFilter(String matchName) {
			if (!matchName.endsWith("\\."))
				matchName += "\\.";
			this.matchName = Pattern.compile(matchName,
					Pattern.CASE_INSENSITIVE);
		}

		public boolean accept(File dir, String name) {
			// Ignore the directory..
			return this.matchName.matcher(name).find();
		}
	}

	private Map<String, Object> input;
	private List<Map<String, Object>> output = new ArrayList<Map<String, Object>>();

	@SuppressWarnings("unchecked")
	public static void main(String[] args) throws YamlException, IOException {
		YamlReader reader = new YamlReader(new FileReader(args[0]));
		Batcher batcher = new Batcher((Map<String, Object>) reader.read());
		batcher.batch();
		batcher.writeOutput();
		// Log
		System.out
				.println("For each Batch one AHRD input-file has been written into: '"
						+ batcher.getInput().get(BATCH_YMLS_DIR_KEY).toString()
						+ "'.");
		System.out
				.println("Created shell-script to start AHRD on all Batches in parallel: '"
						+ batcher.getInput().get(SHELL_SCRIPT_KEY) + "'.");
	}

	public Batcher(Map<String, Object> input) {
		super();
		setInput(input);
	}

	/**
	 * @return The directory concatonated with the found file.
	 */
	public String findFileInDirectory(String directory, String batchName) {
		// Delete file-extension in Batch-Name
		String theBatchName = batchName.replaceAll("\\.\\S+$", "");
		// Handle the / terminating the directory:
		if (!directory.endsWith("/"))
			directory += "/";
		File theDir = new File(directory);
		// Output
		String[] matchingFiles = null;
		if (theDir.isDirectory())
			matchingFiles = theDir.list(new BatchFilenameFilter(theBatchName));
		String matchingFileName = null;
		if (matchingFiles != null && matchingFiles.length == 1)
			matchingFileName = directory + matchingFiles[0];
		else if (matchingFiles != null && matchingFiles.length > 1)
			throw new RuntimeException("Found " + matchingFiles.length
					+ " files in directory '" + directory
					+ "' case-insensitively matching Batch-Name '"
					+ theBatchName + "': "
					+ Arrays.asList(matchingFiles).toString());
		return matchingFileName;
	}

	/**
	 * Argument batchName is expected to be the complete file-name!
	 */
	@SuppressWarnings("unchecked")
	public Map<String, Object> generateYml(String batchName) {
		Map<String, Object> batchYml = new HashMap<String, Object>();

		// Store batchName:
		batchYml.put(BATCH_NAME_KEY, batchName);

		// Set proteins_fasta:
		batchYml.put(Settings.PROTEINS_FASTA_KEY,
				appendSlashIfNotPresent(getInput().get(PROTEINS_DIR_KEY)
						.toString())
						+ batchName);

		// Set references, if given:
		if (getInput().containsKey(REFERENCES_DIR_KEY)) {
			batchYml.put(Settings.REFERENCES_FASTA_KEY, generatePathToFile(
					batchName, REFERENCES_DIR_KEY, null));
		}
		// Set blast2go-results, if given:
		if (getInput().containsKey(BLAST_2_GO_RESULTS_DIR_KEY)) {
			batchYml.put(Settings.BLAST_2_GO_ANNOT_FILE_KEY,
					generatePathToFile(batchName, BLAST_2_GO_RESULTS_DIR_KEY,
							null));
		}
		// Store F-Score-Beta-Parameter, if given:
		if (getInput().containsKey(Settings.F_MEASURE_BETA_PARAM_KEY)) {
			batchYml.put(Settings.F_MEASURE_BETA_PARAM_KEY, getInput().get(
					Settings.F_MEASURE_BETA_PARAM_KEY));
		}

		// Put weight-parameters:
		batchYml.put(Settings.TOKEN_SCORE_BIT_SCORE_WEIGHT, getInput().get(
				Settings.TOKEN_SCORE_BIT_SCORE_WEIGHT));
		batchYml.put(Settings.TOKEN_SCORE_DATABASE_SCORE_WEIGHT, getInput()
				.get(Settings.TOKEN_SCORE_DATABASE_SCORE_WEIGHT));
		batchYml.put(Settings.TOKEN_SCORE_OVERLAP_SCORE_WEIGHT, getInput().get(
				Settings.TOKEN_SCORE_OVERLAP_SCORE_WEIGHT));
		batchYml
				.put(
						Settings.DESCRIPTION_SCORE_RELATIVE_DESCIPTION_FREQUENCY_WEIGHT,
						getInput()
								.get(
										Settings.DESCRIPTION_SCORE_RELATIVE_DESCIPTION_FREQUENCY_WEIGHT));

		// Put data for each Blast-Db into Yml-Hash
		Map<String, Object> blast_dbs = new HashMap<String, Object>();
		Map<String, Object> inputBlastDbs = (Map<String, Object>) getInput()
				.get(Settings.BLAST_DBS_KEY);
		for (String blastDbName : inputBlastDbs.keySet()) {
			// in:
			Map<String, String> inputBlastDb = (Map<String, String>) inputBlastDbs
					.get(blastDbName);
			// out:
			Map<String, String> blastDbYml = new HashMap<String, String>();
			blastDbYml.put(Settings.BLAST_DB_WEIGHT_KEY, inputBlastDb
					.get(Settings.BLAST_DB_WEIGHT_KEY));
			blastDbYml.put(Settings.BLAST_BLACKLIST_KEY, inputBlastDb
					.get(Settings.BLAST_BLACKLIST_KEY));
			blastDbYml.put(Settings.BLAST_FILTER_KEY, inputBlastDb
					.get(Settings.BLAST_FILTER_KEY));
			blastDbYml.put(Settings.TOKEN_BLACKLIST_KEY, inputBlastDb
					.get(Settings.TOKEN_BLACKLIST_KEY));
			// Weight:
			blastDbYml.put(Settings.DESCRIPTION_SCORE_BIT_SCORE_WEIGHT,
					inputBlastDb
							.get(Settings.DESCRIPTION_SCORE_BIT_SCORE_WEIGHT));
			// Find matching Blast-Result-File
			blastDbYml.put(Settings.BLAST_RESULT_FILE_KEY, findFileInDirectory(
					inputBlastDb.get(BLAST_RESULTS_DIR_KEY), batchName));
			blast_dbs.put(blastDbName, blastDbYml);
		}
		batchYml.put(Settings.BLAST_DBS_KEY, blast_dbs);

		// Interpro-Data, only if given:
		if (getInput().get(Settings.INTERPRO_DATABASE_KEY) != null
				&& (getInput().get(INTERPRO_RESULTS_DIR_KEY) != null || getInput()
						.get(INTERPRO_RESULTS_FILE_KEY) != null)) {
			batchYml.put(Settings.INTERPRO_DATABASE_KEY, getInput().get(
					Settings.INTERPRO_DATABASE_KEY));
			batchYml.put(Settings.INTERPRO_RESULT_KEY, generatePathToFile(
					batchName, INTERPRO_RESULTS_DIR_KEY,
					INTERPRO_RESULTS_FILE_KEY));
		}

		// Gene-Ontology-Result, if given:
		if (getInput().get(GENE_ONTOLOGY_RESULTS_DIR_KEY) != null
				|| getInput().get(GENE_ONTOLOGY_RESULTS_FILE_KEY) != null) {
			batchYml.put(Settings.GENE_ONTOLOGY_RESULT_KEY, generatePathToFile(
					batchName, GENE_ONTOLOGY_RESULTS_DIR_KEY,
					GENE_ONTOLOGY_RESULTS_FILE_KEY));
		}

		// Output-File:
		String outputFile = getInput().get(OUTPUT_DIR_KEY).toString();
		if (!outputFile.endsWith("/"))
			outputFile += "/";
		outputFile += batchName.replaceAll("\\.\\S+$", "") + "_ahrd_out.csv";
		batchYml.put(Settings.OUTPUT_KEY, outputFile);

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
			batchYml.put(Settings.TEMPERATURE_KEY, (String) getInput().get(
					Settings.TEMPERATURE_KEY));
		if (getInput().get(Settings.COOL_DOWN_BY_KEY) != null)
			batchYml.put(Settings.COOL_DOWN_BY_KEY, (String) getInput().get(
					Settings.COOL_DOWN_BY_KEY));

		// Path to Batch's YML-File:
		batchYml.put(PATH_TO_BATCH_YML_KEY, generatePathToBatchYml(batchName));
		return batchYml;
	}

	public String generatePathToFile(String batchName, String dirKey,
			String fileKey) {
		String out = null;
		if (getInput().get(dirKey) != null
				&& !((String) getInput().get(dirKey)).equals("")) {
			out = findFileInDirectory(getInput().get(dirKey).toString(),
					batchName);
		} else if (getInput().get(fileKey) != null
				&& !((String) getInput().get(fileKey)).equals("")) {
			out = (String) getInput().get(fileKey);
		} else {
			throw new IllegalArgumentException("Missing either " + dirKey
					+ " or " + fileKey + " in batcher_input.yml");
		}
		return out;
	}

	public String generateAhrdCall(String pathToBatchYml, String batchName) {
		String ahrdCall = getInput().get(AHRD_CALL_KEY).toString().replaceAll(
				AHRD_CALL_BATCH, pathToBatchYml);
		return ahrdCall.replaceAll(AHRD_CALL_BATCH_NAME, batchName);
	}

	public String generatePathToBatchYml(String batchName) {
		String pathToBatchYml = getInput().get(BATCH_YMLS_DIR_KEY).toString();
		if (!pathToBatchYml.endsWith("/"))
			pathToBatchYml += "/";
		pathToBatchYml += batchName.replaceAll("\\.\\S+$", ".yml");
		return pathToBatchYml;
	}

	public void batch() {
		File proteinsDir = new File(getInput().get(PROTEINS_DIR_KEY).toString());
		if (proteinsDir.isDirectory()) {
			for (File proteinFile : proteinsDir.listFiles()) {
				if (!proteinFile.isDirectory()) {
					getOutput().add(generateYml(proteinFile.getName()));
				}
			}
		}
	}

	public void writeOutput() throws IOException {
		String batchYmlsDir = getInput().get(BATCH_YMLS_DIR_KEY).toString();
		if (!batchYmlsDir.endsWith("/"))
			batchYmlsDir += "/";

		BufferedWriter shellScriptBw = new BufferedWriter(new FileWriter(
				getInput().get(SHELL_SCRIPT_KEY).toString()));

		for (Map<String, Object> batchYml : getOutput()) {
			// delete the reference to the Batch's name, as we do not want to
			// see it in the output-yml:
			String batchName = batchYml.remove(BATCH_NAME_KEY).toString();
			// delete the reference in Batch's YML to itself
			String pathToBatchYml = batchYml.remove(PATH_TO_BATCH_YML_KEY)
					.toString();
			YamlWriter writer = new YamlWriter(new FileWriter(pathToBatchYml));
			writer.write(batchYml);
			writer.close();
			// generate and store Shell-Script-Commands:
			shellScriptBw.write(generateAhrdCall(pathToBatchYml, batchName)
					+ "\n");
		}

		// Write Shell-Script:
		shellScriptBw.close();
	}

	public String appendSlashIfNotPresent(String inDirPath) {
		if (!inDirPath.endsWith("/"))
			inDirPath += "/";
		return inDirPath;
	}

	/**
	 * Gets the input for this instance.
	 * 
	 * @return The input.
	 */
	public Map<String, Object> getInput() {
		return this.input;
	}

	/**
	 * Sets the input for this instance.
	 * 
	 * @param input
	 *            The input.
	 */
	public void setInput(Map<String, Object> input) {
		this.input = input;
	}

	/**
	 * Gets the output for this instance.
	 * 
	 * @return The output.
	 */
	public List<Map<String, Object>> getOutput() {
		return this.output;
	}

	/**
	 * Sets the output for this instance.
	 * 
	 * @param output
	 *            The output.
	 */
	public void setOutput(List<Map<String, Object>> output) {
		this.output = output;
	}

}
