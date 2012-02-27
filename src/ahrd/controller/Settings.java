package ahrd.controller;

import static ahrd.controller.Utils.fromFile;
import static ahrd.controller.Utils.readFile;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.esotericsoftware.yamlbeans.YamlReader;

/**
 * Using the <em>thread-local</em> singleton-pattern to hold in a central place
 * all AHRD's current settings. This eases greatly access of input-values from
 * everywhere.
 * 
 * @author hallab, klee
 */
public class Settings implements Cloneable {

	/**
	 * Thread-Local Singleton of the current AHRD-Run's settings:
	 */
	private static final ThreadLocal<Settings> settings = new ThreadLocal<Settings>();

	public static Settings getSettings() {
		return settings.get();
	}

	public static void setSettings(Settings s) {
		settings.set(s);
	}

	/**
	 * Keys to parse the YML-Input:
	 */
	public static final String PROTEINS_FASTA_KEY = "proteins_fasta";
	public static final String BLAST_DBS_KEY = "blast_dbs";
	public static final String BLAST_DB_WEIGHT_KEY = "weight";
	public static final String BLAST_RESULT_FILE_KEY = "file";
	public static final String BLAST_BLACKLIST_KEY = "blacklist";
	public static final String BLAST_FILTER_KEY = "filter";
	public static final String TOKEN_BLACKLIST_KEY = "token_blacklist";
	public static final String INTERPRO_DATABASE_KEY = "interpro_database";
	public static final String INTERPRO_RESULT_KEY = "interpro_result";
	public static final String GENE_ONTOLOGY_RESULT_KEY = "gene_ontology_result";
	public static final String OUTPUT_KEY = "output";
	public static final String WRITE_SCORES_TO_OUTPUT = "write_scores_to_output";
	public static final String WRITE_BEST_BLAST_HITS_TO_OUTPUT = "write_best_blast_hits_to_output";
	public static final String WRITE_TOKEN_SET_TO_OUTPUT = "write_token_set_to_output";
	public static final String TOKEN_SCORE_BIT_SCORE_WEIGHT = "token_score_bit_score_weight";
	public static final String TOKEN_SCORE_DATABASE_SCORE_WEIGHT = "token_score_database_score_weight";
	public static final String TOKEN_SCORE_OVERLAP_SCORE_WEIGHT = "token_score_overlap_score_weight";
	public static final String DESCRIPTION_SCORE_RELATIVE_DESCIPTION_FREQUENCY_WEIGHT = "description_score_relative_description_frequency_weight";
	public static final String DESCRIPTION_SCORE_BIT_SCORE_WEIGHT = "description_score_bit_score_weight";
	public static final String REFERENCES_FASTA_KEY = "references_fasta";
	public static final String F_MEASURE_BETA_PARAM_KEY = "f_measure_beta_parameter";
	public static final String BLAST_2_GO_ANNOT_FILE_KEY = "blast2go";
	public static final Double PERCENTAGE_MUTATOR_SEED = 0.10;
	public static final Integer BLAST_DB_WEIGHT_MUTATOR_SEED = 10;
	public static final String TEMPERATURE_KEY = "temperature";
	public static final String COOL_DOWN_BY_KEY = "cool_down_by";
	public static final String NO_START_POSITIONS_IN_PARAM_SPACE = "no_start_positions_in_parameter_space";
	public static final String REMEMBER_SIMULATED_ANNEALING_PATH_KEY = "remember_simulated_annealing_path";

	/**
	 * Fields:
	 */
	private String pathToProteinsFasta;
	private String pathToReferencesFasta;
	private String pathToInterproDatabase;
	private String pathToInterproResults;
	private String pathToGeneOntologyResults;
	private String pathToOutput;
	/**
	 * Parameters representing weights and factors in the various formulas used
	 * in AHRD. They are subject to optimization and can be set by the user.
	 */
	private Parameters parameters = new Parameters();
	private Boolean writeTokenSetToOutput;
	private Boolean writeBestBlastHitsToOutput;
	private Boolean writeScoresToOutput;
	/**
	 * F-Measure's Beta-Parameter as set in the input.yml or default 1.0
	 */
	private Double fMeasureBetaParameter = 1.0;
	private Map<String, Map<String, String>> blastDbSettings = new HashMap<String, Map<String, String>>();
	private Map<String, List<String>> blastResultsBlacklists = new HashMap<String, List<String>>();
	private Map<String, List<String>> blastResultsFilter = new HashMap<String, List<String>>();
	private Map<String, List<String>> tokenBlacklists = new HashMap<String, List<String>>();
	private String pathToBlast2GoAnnotations;
	/**
	 * For the <strong>simulated annealing</strong> algorithm, this will be
	 * current temperature. (Default is 1000)
	 */
	private Integer temperature = 1000;
	/**
	 * For the <strong>simulated annealing</strong> algorithm, this will be
	 * value the current temperature gets cooled down each step. (Default is 1)
	 */
	private Integer coolDownBy = 1;
	/**
	 * Break with the classic simulated annealing approach and remember each
	 * visited Parameter-Set and its score. This enables speeding up the
	 * optimization with the drawback of higher memory usage.
	 */
	private boolean rememberSimulatedAnnealingPath = false;

	/**
	 * Construct from contents of file 'AHRD_input.yml'.
	 * 
	 * @throws IOException
	 */
	public Settings(String pathToYml) throws IOException {
		super();
		this.initialize(pathToYml);
	}

	/**
	 * Initializes an Instance with content read from a YML-File:
	 * 
	 * @throws IOException
	 */
	@SuppressWarnings("unchecked")
	public void initialize(String pathToYml) throws IOException {
		YamlReader reader = new YamlReader(new FileReader(pathToYml));
		Map<String, Object> input = (Map<String, Object>) reader.read();
		this.blastDbSettings = (Map<String, Map<String, String>>) input
				.get(BLAST_DBS_KEY);
		setPathToProteinsFasta((String) input.get(PROTEINS_FASTA_KEY));
		setPathToInterproDatabase((String) input.get(INTERPRO_DATABASE_KEY));
		setPathToInterproResults((String) input.get(INTERPRO_RESULT_KEY));
		setPathToGeneOntologyResults((String) input
				.get(GENE_ONTOLOGY_RESULT_KEY));
		setPathToOutput((String) input.get(OUTPUT_KEY));
		setTokenScoreBitScoreWeight(Double.parseDouble((String) input
				.get(TOKEN_SCORE_BIT_SCORE_WEIGHT)));
		setTokenScoreDatabaseScoreWeight(Double.parseDouble((String) input
				.get(TOKEN_SCORE_DATABASE_SCORE_WEIGHT)));
		setTokenScoreOverlapScoreWeight(Double.parseDouble((String) input
				.get(TOKEN_SCORE_OVERLAP_SCORE_WEIGHT)));
		setDescriptionScorePatternFactorWeight(Double
				.parseDouble((String) input
						.get(DESCRIPTION_SCORE_RELATIVE_DESCIPTION_FREQUENCY_WEIGHT)));
		setWriteTokenSetToOutput(Boolean.parseBoolean((String) input
				.get(WRITE_TOKEN_SET_TO_OUTPUT)));
		setWriteBestBlastHitsToOutput(Boolean.parseBoolean((String) input
				.get(WRITE_BEST_BLAST_HITS_TO_OUTPUT)));
		setWriteScoresToOutput(Boolean.parseBoolean((String) input
				.get(WRITE_SCORES_TO_OUTPUT)));
		// Generate the Blacklists and Filters for each Blast-Database from
		// their appropriate files:
		for (String blastDatabaseName : getBlastDatabases()) {
			this.blastResultsBlacklists
					.put(blastDatabaseName,
							fromFile(getPathToBlastResultsBlackList(blastDatabaseName)));
			this.blastResultsFilter.put(blastDatabaseName,
					fromFile(getPathToBlastResultsFilter(blastDatabaseName)));
			this.tokenBlacklists.put(blastDatabaseName,
					fromFile(getPathToTokenBlacklist(blastDatabaseName)));
			// Set Database-Weights and Description-Score-Bit-Score-Weight:
			this.getParameters().setBlastDbWeight(
					blastDatabaseName,
					this.getBlastDbSettings(blastDatabaseName).get(
							Settings.BLAST_DB_WEIGHT_KEY));
			this.getParameters().setDescriptionScoreBitScoreWeight(
					blastDatabaseName,
					this.getBlastDbSettings(blastDatabaseName).get(
							Settings.DESCRIPTION_SCORE_BIT_SCORE_WEIGHT));
		}
		// If started to train the algorithm references are stored in this file:
		setPathToReferencesFasta((String) input.get(REFERENCES_FASTA_KEY));
		// If started in training-mode the F-Measure's Beta-Parameter can be set
		// to some other value than 1.0
		if (input.get(F_MEASURE_BETA_PARAM_KEY) != null)
			this.fMeasureBetaParameter = Double.parseDouble((String) input
					.get(F_MEASURE_BETA_PARAM_KEY));
		// If started to compare AHRD with Blast2Go, enable reading of
		// B2G-Annotations:
		if (input.get(BLAST_2_GO_ANNOT_FILE_KEY) != null)
			this.pathToBlast2GoAnnotations = (String) input
					.get(BLAST_2_GO_ANNOT_FILE_KEY);
		// Simulated Annealing can be started with custom temperature and value
		// it is cooled-down by each step:
		if (input.get(TEMPERATURE_KEY) != null)
			setTemperature(Integer
					.parseInt((String) input.get(TEMPERATURE_KEY)));
		if (input.get(COOL_DOWN_BY_KEY) != null)
			this.coolDownBy = Integer.parseInt((String) input
					.get(COOL_DOWN_BY_KEY));
		if (input.get(REMEMBER_SIMULATED_ANNEALING_PATH_KEY) != null
				&& Boolean.parseBoolean(input.get(
						REMEMBER_SIMULATED_ANNEALING_PATH_KEY).toString()))
			this.rememberSimulatedAnnealingPath = true;
	}

	/**
	 * Returns a clone of this instance. <strong>Only</strong> all primitive
	 * fields and the Blast-Database-Parameters are actually cloned. All other
	 * fields still refer to <strong>the same objects</strong>.
	 * <em>So be very careful using this method.</em> It has been written in
	 * this manner to fulfill requirements and minimize memory-usage.
	 */
	public Settings clone() {
		Settings clone;
		try {
			clone = (Settings) super.clone();
		} catch (CloneNotSupportedException e) {
			e.printStackTrace(System.err);
			return null;
		}
		// Clone the Blast-Database-Settings-Map and it's Parameters:
		Map<String, Map<String, String>> blastDbSettings = new HashMap<String, Map<String, String>>();
		for (String blastDb : getBlastDbSettings().keySet()) {
			blastDbSettings.put(blastDb, new HashMap<String, String>());
			for (String iterKey : getBlastDbSettings(blastDb).keySet()) {
				blastDbSettings.get(blastDb).put(new String(iterKey),
						new String(getBlastDbSettings(blastDb).get(iterKey)));
			}
		}
		clone.setBlastDbSettings(blastDbSettings);
		// Clone the Parameters subject to optimization:
		clone.setParameters(this.getParameters().clone());
		return clone;
	}

	public boolean hasValidInterproDatabaseAndResultFile() {
		if (getPathToInterproDatabase() == null
				|| getPathToInterproResults() == null)
			return false;
		// ELSE:
		File iprDb = new File(getPathToInterproDatabase());
		File iprRes = new File(getPathToInterproResults());
		return (iprDb.canRead() && iprDb.length() > 0 && iprRes.canRead() && iprRes
				.length() > 0);
	}

	/**
	 * Break with the classic simulated annealing approach and remember each
	 * visited Parameter-Set and its score. This enables speeding up the
	 * optimization with the drawback of higher memory usage.
	 * 
	 * @return boolean - flag
	 */
	public boolean rememberSimulatedAnnealingPath() {
		return this.rememberSimulatedAnnealingPath;
	}

	/**
	 * @param blastDatabaseName
	 * @return Map<String, String> the Settings for the argument Blast-Database
	 */
	private Map<String, String> getBlastDbSettings(String blastDatabaseName) {
		return this.blastDbSettings.get(blastDatabaseName);
	}

	/**
	 * @return Set<String> the names of the blast-databases used in the current
	 *         AHRD-Run.
	 */
	public Set<String> getBlastDatabases() {
		return this.blastDbSettings.keySet();
	}

	public Integer getBlastDbWeight(String blastDatabaseName) {
		return getParameters().getBlastDbWeight(blastDatabaseName);
	}

	public void setBlastDbWeight(String blastDatabaseName, String bdbw) {
		getParameters().setBlastDbWeight(blastDatabaseName, bdbw);
	}

	public Double getDescriptionScoreBitScoreWeight(String blastDatabaseName) {
		return getParameters().getDescriptionScoreBitScoreWeight(
				blastDatabaseName);
	}

	public void setDescriptionScoreBitScoreWeight(String blastDatabaseName,
			String dsbsw) {
		getParameters().setDescriptionScoreBitScoreWeight(blastDatabaseName,
				dsbsw);
	}

	public String getPathToBlastResults(String blastDatabaseName) {
		return getBlastDbSettings(blastDatabaseName).get(BLAST_RESULT_FILE_KEY);
	}

	private String getPathToBlastResultsBlackList(String blastDatabaseName) {
		return getBlastDbSettings(blastDatabaseName).get(BLAST_BLACKLIST_KEY);
	}

	public List<String> getBlastResultsBlackList(String blastDatabaseName) {
		return this.blastResultsBlacklists.get(blastDatabaseName);
	}

	private String getPathToBlastResultsFilter(String blastDatabaseName) {
		return getBlastDbSettings(blastDatabaseName).get(BLAST_FILTER_KEY);
	}

	public List<String> getBlastResultsFilter(String blastDatabaseName) {
		return this.blastResultsFilter.get(blastDatabaseName);
	}

	private String getPathToTokenBlacklist(String blastDatabaseName) {
		return getBlastDbSettings(blastDatabaseName).get(TOKEN_BLACKLIST_KEY);
	}

	public List<String> getTokenBlackList(String blastDatabaseName) {
		return this.tokenBlacklists.get(blastDatabaseName);
	}

	private String getPathToProteinsFasta() {
		return pathToProteinsFasta;
	}

	public String getProteinsFasta() throws IOException {
		return readFile(getPathToProteinsFasta());
	}

	private String getPathToReferencesFasta() {
		return pathToReferencesFasta;
	}

	public String getReferencesFasta() throws IOException {
		return readFile(getPathToReferencesFasta());
	}

	public void setPathToProteinsFasta(String pathToProteinsFasta) {
		this.pathToProteinsFasta = pathToProteinsFasta;
	}

	public String getPathToInterproDatabase() {
		return pathToInterproDatabase;
	}

	public void setPathToInterproDatabase(String pathToInterproDatabase) {
		this.pathToInterproDatabase = pathToInterproDatabase;
	}

	public String getPathToInterproResults() {
		return pathToInterproResults;
	}

	public boolean hasInterproAnnotations() {
		return getPathToInterproDatabase() != null
				&& (new File(getPathToInterproDatabase())).exists()
				&& getPathToInterproResults() != null
				&& (new File(getPathToInterproResults())).exists();
	}

	public void setPathToInterproResults(String pathToInterproResults) {
		this.pathToInterproResults = pathToInterproResults;
	}

	public String getPathToGeneOntologyResults() {
		return pathToGeneOntologyResults;
	}

	public boolean hasGeneOntologyAnnotations() {
		return getPathToGeneOntologyResults() != null
				&& (new File(getPathToGeneOntologyResults())).exists();
	}

	public void setPathToGeneOntologyResults(String pathToGeneOntologyResults) {
		this.pathToGeneOntologyResults = pathToGeneOntologyResults;
	}

	public String getPathToOutput() {
		return pathToOutput;
	}

	public void setPathToOutput(String pathToOutput) {
		this.pathToOutput = pathToOutput;
	}

	public Double getTokenScoreBitScoreWeight() {
		return getParameters().getTokenScoreBitScoreWeight();
	}

	public void setTokenScoreBitScoreWeight(Double tokenScoreBitScoreWeight) {
		this.getParameters().setTokenScoreBitScoreWeight(
				tokenScoreBitScoreWeight);
	}

	public Double getTokenScoreDatabaseScoreWeight() {
		return getParameters().getTokenScoreDatabaseScoreWeight();
	}

	public void setTokenScoreDatabaseScoreWeight(
			Double tokenScoreDatabaseScoreWeight) {
		this.getParameters().setTokenScoreDatabaseScoreWeight(
				tokenScoreDatabaseScoreWeight);
	}

	public Double getTokenScoreOverlapScoreWeight() {
		return getParameters().getTokenScoreOverlapScoreWeight();
	}

	public void setTokenScoreOverlapScoreWeight(
			Double tokenScoreOverlapScoreWeight) {
		this.getParameters().setTokenScoreOverlapScoreWeight(
				tokenScoreOverlapScoreWeight);
	}

	public Double getDescriptionScorePatternFactorWeight() {
		return getParameters().getDescriptionScorePatternFactorWeight();
	}

	public void setDescriptionScorePatternFactorWeight(
			Double descriptionScorePatternFactorWeight) {
		this.getParameters().setDescriptionScorePatternFactorWeight(
				descriptionScorePatternFactorWeight);
	}

	public Boolean getWriteTokenSetToOutput() {
		return writeTokenSetToOutput;
	}

	public void setWriteTokenSetToOutput(Boolean writeTokenSetToOutput) {
		this.writeTokenSetToOutput = writeTokenSetToOutput;
	}

	public Boolean getWriteBestBlastHitsToOutput() {
		return writeBestBlastHitsToOutput;
	}

	public void setWriteBestBlastHitsToOutput(Boolean writeBestBlastHitsToOutput) {
		this.writeBestBlastHitsToOutput = writeBestBlastHitsToOutput;
	}

	public Boolean getWriteScoresToOutput() {
		return writeScoresToOutput;
	}

	public void setWriteScoresToOutput(Boolean writeScoresToOutput) {
		this.writeScoresToOutput = writeScoresToOutput;
	}

	public Map<String, Map<String, String>> getBlastDbSettings() {
		return blastDbSettings;
	}

	public void setBlastDbSettings(
			Map<String, Map<String, String>> blastDbSettings) {
		this.blastDbSettings = blastDbSettings;
	}

	public void setPathToReferencesFasta(String pathToReferencesFasta) {
		this.pathToReferencesFasta = pathToReferencesFasta;
	}

	public boolean isInTrainingMode() {
		return (getPathToReferencesFasta() != null && getPathToReferencesFasta() != "");
	}

	/**
	 * @return Double - F-Measure's Beta-Parameter as set in the input.yml or
	 *         default 1.0
	 */
	public Double getFMeasureBetaParameter() {
		return this.fMeasureBetaParameter;
	}

	public String getPathToBlast2GoAnnotations() {
		return pathToBlast2GoAnnotations;
	}

	public List<String> getBlast2GoAnnotations() throws IOException {
		return fromFile(getPathToBlast2GoAnnotations());
	}

	public Double getAvgEvaluationScore() {
		return getParameters().getAvgEvaluationScore();
	}

	public void setAvgEvaluationScore(Double avgEvaluationScore) {
		this.getParameters().setAvgEvaluationScore(avgEvaluationScore);
	}

	public Double getAvgTruePositivesRate() {
		return getParameters().getAvgTruePositivesRate();
	}

	public void setAvgTruePositivesRate(Double avgTruePositivesRate) {
		this.getParameters().setAvgTruePositivesRate(avgTruePositivesRate);
	}

	public Double getAvgFalsePositivesRate() {
		return getParameters().getAvgFalsePositivesRate();
	}

	public void setAvgFalsePositivesRate(Double avgFalsePositivesRate) {
		this.getParameters().setAvgFalsePositivesRate(avgFalsePositivesRate);
	}

	public Integer getTemperature() {
		return temperature;
	}

	public void setTemperature(Integer temperature) {
		this.temperature = temperature;
	}

	public Integer getCoolDownBy() {
		return coolDownBy;
	}

	public Parameters getParameters() {
		return parameters;
	}

	public void setParameters(Parameters parameters) {
		this.parameters = parameters;
	}
}
