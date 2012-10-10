package ahrd.controller;

import static ahrd.controller.Utils.fromFile;
import static ahrd.controller.Utils.readFile;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
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
	public static final String DOMAIN_WEIGHTS_DATABASE = "domain_weights_database";
	public static final String DOMAIN_WEIGHTS_POSITION_KEY = "domain_weights_table_position";
	public static final String COMPUTE_DOMAIN_SIMILARITY_ON_KEY = "compute_domain_similarity_on";
	public static final String GENE_ONTOLOGY_RESULT_KEY = "gene_ontology_result";
	public static final String OUTPUT_KEY = "output";
	public static final String SIMULATED_ANNEALING_PATH_LOG_KEY = "path_log";
	public static final String WRITE_SCORES_TO_OUTPUT = "write_scores_to_output";
	public static final String WRITE_BEST_BLAST_HITS_TO_OUTPUT = "write_best_blast_hits_to_output";
	public static final String WRITE_TOKEN_SET_TO_OUTPUT = "write_token_set_to_output";
	public static final String HRD_SCORES_OUTPUT_PATH = "hrd_scores_output";
	public static final String TOKEN_SCORE_BIT_SCORE_WEIGHT = "token_score_bit_score_weight";
	public static final String TOKEN_SCORE_DATABASE_SCORE_WEIGHT = "token_score_database_score_weight";
	public static final String TOKEN_SCORE_OVERLAP_SCORE_WEIGHT = "token_score_overlap_score_weight";
	public static final String DESCRIPTION_SCORE_RELATIVE_DESCIPTION_FREQUENCY_WEIGHT = "description_score_relative_description_frequency_weight";
	public static final String DESCRIPTION_SCORE_BIT_SCORE_WEIGHT = "description_score_bit_score_weight";
	public static final String REFERENCES_FASTA_KEY = "references_fasta";
	public static final String F_MEASURE_BETA_PARAM_KEY = "f_measure_beta_parameter";
	public static final String BLAST_2_GO_ANNOT_FILE_KEY = "blast2go";
	public static final String TEMPERATURE_KEY = "temperature";
	public static final String COOL_DOWN_BY_KEY = "cool_down_by";
	public static final String OPTIMIZATION_ACCEPTANCE_PROBABILITY_SCALING_FACTOR_KEY = "optimization_acceptance_probability_scaling_factor";
	public static final String MUTATOR_MEAN_KEY = "mutator_mean";
	public static final String MUTATOR_DEVIATION_KEY = "mutator_deviation";
	public static final String NO_START_POSITIONS_IN_PARAM_SPACE = "no_start_positions_in_parameter_space";
	public static final String REMEMBER_SIMULATED_ANNEALING_PATH_KEY = "remember_simulated_annealing_path";
	public static final String P_MUTATE_SAME_PARAMETER_SCALE_KEY = "p_mutate_same_parameter_scale";
	public static final String FIND_HIGHEST_POSSIBLE_EVALUATION_SCORE_KEY = "find_highest_possible_evaluation_score";
	public static final String OUTPUT_FASTA_KEY = "output_fasta";
	public static final String DESCRIPTION_SCORE_DOMAIN_SIMILARITY_WEIGHT_KEY = "description_score_domain_similarity_weight";
	public static final String TOKEN_SCORE_DOMAIN_SIMILARITY_WEIGHT_KEY = "token_score_domain_similarity_weight";
	public static final String WRITE_DOMAIN_ARCHITECTURE_SIMILARITY_SCORES_TO_OUTPUT = "write_domain_architecture_similarity_scores_to_output";

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
	 * File to write the AHRD-Scores of each BlastHit's Description into, if
	 * requested.
	 */
	private String pathToHRDScoresOutput;
	/**
	 * Trainer logs path through parameter- and score-space into this file:
	 */
	private String pathToSimulatedAnnealingPathLog;
	/**
	 * Parameters representing weights and factors in the various formulas used
	 * in AHRD. They are subject to optimization and can be set by the user.
	 */
	private Parameters parameters = new Parameters();
	private Boolean writeTokenSetToOutput;
	private Boolean writeBestBlastHitsToOutput;
	/**
	 * Forces AHRD to write out all internal scores (Sum(Token-Scores),
	 * Description- and Lexical-Scores, etc.
	 */
	private Boolean writeScoresToOutput;
	/**
	 * F-Measure's Beta-Parameter as set in the input.yml or default 1.0
	 */
	private Double fMeasureBetaParameter = 1.0;
	private Map<String, Map<String, String>> blastDbSettings = new HashMap<String, Map<String, String>>();
	private List<String> sortedBlastDatabaseNames;
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
	 * For the <strong>simulated annealing</strong> algorithm, this is the
	 * scaling factor for the probability distribution P('Accept worse scoring
	 * Parameter-Set') := exp(-delta_scores*scaling_factor/current-temperature).
	 */
	private Double optimizationAcceptanceProbabilityScalingFactor = 200000000.0;
	/**
	 * In simulated annealing optimization each cycle has to mutate the current
	 * Parameter-Set to generate a neighboring set in parameter space. The
	 * random value used to add or subtract to a single parameter is Gaussian
	 * distributed and has the following mean:
	 */
	private Double mutatorMean = 0.2;
	/**
	 * In simulated annealing optimization each cycle has to mutate the current
	 * Parameter-Set to generate a neighboring set in parameter space. The
	 * random value used to add or subtract to a single parameter is Gaussian
	 * distributed and has the following standard deviation:
	 */
	private Double mutatorDeviation = 0.25;
	/**
	 * If the last optimization step was done with better performing parameters,
	 * randomly decide to mutate the same Parameter to generate a new Neighbor
	 * in Parameter-Space. By this simulated annealing walks more likely uphill
	 * in Parameter-Score-Space. The probability P('Mutate same Parameter') :=
	 * 0, if score was not increased, (exp(1-increase.score)+s)/(exp(0)+s) else
	 * 
	 * This is the scaling parameter s in above formula.
	 */
	private Double pMutateSameParameterScale = 0.7;
	/**
	 * Break with the classic simulated annealing approach and remember each
	 * visited Parameter-Set and its score. This enables speeding up the
	 * optimization with the drawback of higher memory usage.
	 */
	private boolean rememberSimulatedAnnealingPath = false;
	/**
	 * Evaluation or Optimization might be interested in the highest possibly
	 * achievable evaluation-score:
	 */
	private boolean findHighestPossibleEvaluationScore = false;
	/**
	 * Write output as fasta-file?
	 */
	private boolean outputFasta = false;
	/**
	 * Path to Domain-Weight file as downloadable from
	 * http://pat.kobic.re.kr/wdac/data/domain_scores.gz
	 */
	private String pathToDomainWeightsDatabase;
	/**
	 * Path to interproscan results of the Proteins referenced in the BlastHits.
	 * This data is needed for the Domain-Scoring.
	 */
	private String pathToInterproResults4BlastHits;
	/**
	 * Path to Pfam results of the Proteins ireferenced in BlastHts. This data
	 * is needed for the Domain-Scoring.
	 */
	private String pathToPfamResults4BlastHits;

	/**
	 * The configurable weight for the fraction a BlastResult's domain weight
	 * similarity score is going to assume in the final token score.
	 */
	private Double tokenScoreDomainSimilarityWeight = 0.0;
	/**
	 * The configurable weight for the fraction a BlastResult's domain weight
	 * similarity score is going to assume in the final description score.
	 */
	private Double descriptionScoreDomainSimilarityWeight = 0.0;
	/**
	 * If AHRD is to consider similarity of domain architectures it needs to
	 * know, if to base this scoring on annotated Pfam or InterPro domains.
	 */
	private String computeDomainSimilarityOn = null;
	/**
	 * AHRD run in domain_architecture_similarity mode can be asked to append
	 * those scores to the output.
	 */
	private boolean writeDomainArchitectureSimilarityScoresToOutput = false;
	/**
	 * Position of tab delimited table to look up the domain weight:
	 */
	private Integer domainWeightTablePosition = 7;

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
		setPathToDomainWeightsDatabase((String) input
				.get(DOMAIN_WEIGHTS_DATABASE));
		if (input.get(DOMAIN_WEIGHTS_POSITION_KEY) != null)
			setDomainWeightTablePosition(Integer.parseInt((String) input
					.get(DOMAIN_WEIGHTS_POSITION_KEY)));
		setComputeDomainSimilarityOn((String) input
				.get(COMPUTE_DOMAIN_SIMILARITY_ON_KEY));
		if (input.get(TOKEN_SCORE_DOMAIN_SIMILARITY_WEIGHT_KEY) != null)
			setTokenScoreDomainSimilarityWeight(Double
					.parseDouble((String) input
							.get(TOKEN_SCORE_DOMAIN_SIMILARITY_WEIGHT_KEY)));
		if (input.get(DESCRIPTION_SCORE_DOMAIN_SIMILARITY_WEIGHT_KEY) != null)
			setDescriptionScoreDomainSimilarityWeight(Double
					.parseDouble((String) input
							.get(DESCRIPTION_SCORE_DOMAIN_SIMILARITY_WEIGHT_KEY)));
		setPathToGeneOntologyResults((String) input
				.get(GENE_ONTOLOGY_RESULT_KEY));
		setPathToOutput((String) input.get(OUTPUT_KEY));
		if (input.get(HRD_SCORES_OUTPUT_PATH) != null
				&& !input.get(HRD_SCORES_OUTPUT_PATH).equals(""))
			setPathToHRDScoresOutput((String) input.get(HRD_SCORES_OUTPUT_PATH));
		// Trainer logs path through parameter-space here:
		if (input.get(SIMULATED_ANNEALING_PATH_LOG_KEY) != null)
			setPathToSimulatedAnnealingPathLog((String) input
					.get(SIMULATED_ANNEALING_PATH_LOG_KEY));
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
		setWriteDomainArchitectureSimilarityScoresToOutput(Boolean
				.parseBoolean((String) input
						.get(WRITE_DOMAIN_ARCHITECTURE_SIMILARITY_SCORES_TO_OUTPUT)));
		setOutputFasta(Boolean.parseBoolean((String) input
				.get(OUTPUT_FASTA_KEY)));
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
		if (input.get(OPTIMIZATION_ACCEPTANCE_PROBABILITY_SCALING_FACTOR_KEY) != null)
			setOptimizationAcceptanceProbabilityScalingFactor(Double
					.parseDouble((String) input
							.get(OPTIMIZATION_ACCEPTANCE_PROBABILITY_SCALING_FACTOR_KEY)));
		if (input.get(MUTATOR_MEAN_KEY) != null)
			setMutatorMean(Double.parseDouble((String) input
					.get(MUTATOR_MEAN_KEY)));
		if (input.get(MUTATOR_DEVIATION_KEY) != null)
			setMutatorDeviation(Double.parseDouble((String) input
					.get(MUTATOR_DEVIATION_KEY)));
		if (input.get(REMEMBER_SIMULATED_ANNEALING_PATH_KEY) != null
				&& Boolean.parseBoolean(input.get(
						REMEMBER_SIMULATED_ANNEALING_PATH_KEY).toString()))
			this.rememberSimulatedAnnealingPath = true;
		if (input.get(P_MUTATE_SAME_PARAMETER_SCALE_KEY) != null)
			setpMutateSameParameterScale(Double.parseDouble((String) input
					.get(P_MUTATE_SAME_PARAMETER_SCALE_KEY)));
		// Evaluation or Optimization might be interested in the highest
		// possibly achievable evaluation-score:
		if (input.get(FIND_HIGHEST_POSSIBLE_EVALUATION_SCORE_KEY) != null
				&& Boolean.parseBoolean(input.get(
						FIND_HIGHEST_POSSIBLE_EVALUATION_SCORE_KEY).toString()))
			this.findHighestPossibleEvaluationScore = true;
		// Validate input parameters:
		validateComputeDomainSimilarityOn();
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
	 * Validates field <em>computeDomainSimilarityOn</em> is one of the
	 * following:
	 * <ul>
	 * <li>NULL</li>
	 * <li>interpro</li>
	 * <li>pfam</li>
	 * </ul>
	 * 
	 * @throws IllegalArgumentException
	 *             if field is invalid.
	 */
	public void validateComputeDomainSimilarityOn() {
		if (getComputeDomainSimilarityOn() != null
				&& !(getComputeDomainSimilarityOn()
						.equalsIgnoreCase("interpro") || getComputeDomainSimilarityOn()
						.equalsIgnoreCase("pfam"))) {
			throw new IllegalArgumentException(
					"Input parameter '"
							+ COMPUTE_DOMAIN_SIMILARITY_ON_KEY
							+ "' has to be one of the following: NULL, 'interpro' or 'pfam', but was: "
							+ getComputeDomainSimilarityOn());
		}
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

	/**
	 * @return List<String> - The alphabetically sorted List of
	 *         Blast-Database-Names.
	 */
	public List<String> getSortedBlastDatabases() {
		// Only sort ONCE:
		if (this.sortedBlastDatabaseNames == null) {
			this.sortedBlastDatabaseNames = new ArrayList<String>(
					getBlastDatabases());
			Collections.sort(this.sortedBlastDatabaseNames);
		}
		return this.sortedBlastDatabaseNames;
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
	 * Only compute domain similarity scores, if and only if all required input
	 * data is present.
	 * 
	 * @return boolean
	 */
	public boolean isToComputeDomainSimilarities() {
		return (hasInterproAnnotations()
				&& getPathToDomainWeightsDatabase() != null && !getPathToDomainWeightsDatabase()
				.equals(""));
	}

	/**
	 * Evaluation or Optimization might be interested in the highest possibly
	 * achievable evaluation-score.
	 */
	public boolean doFindHighestPossibleEvaluationScore() {
		return this.findHighestPossibleEvaluationScore;
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

	public boolean doOutputFasta() {
		return outputFasta;
	}

	public void setOutputFasta(boolean outputFasta) {
		this.outputFasta = outputFasta;
	}

	public Double getOptimizationAcceptanceProbabilityScalingFactor() {
		return optimizationAcceptanceProbabilityScalingFactor;
	}

	public Double getMutatorMean() {
		return mutatorMean;
	}

	public Double getMutatorDeviation() {
		return mutatorDeviation;
	}

	public void setMutatorMean(Double mutatorMean) {
		this.mutatorMean = mutatorMean;
	}

	public void setMutatorDeviation(Double mutatorDeviation) {
		this.mutatorDeviation = mutatorDeviation;
	}

	public void setOptimizationAcceptanceProbabilityScalingFactor(
			Double optimizationAcceptanceProbabilityScalingFactor) {
		this.optimizationAcceptanceProbabilityScalingFactor = optimizationAcceptanceProbabilityScalingFactor;
	}

	public String getPathToSimulatedAnnealingPathLog() {
		return pathToSimulatedAnnealingPathLog;
	}

	public void setPathToSimulatedAnnealingPathLog(
			String pathToSimulatedAnnealingPathLog) {
		this.pathToSimulatedAnnealingPathLog = pathToSimulatedAnnealingPathLog;
	}

	public Double getpMutateSameParameterScale() {
		return pMutateSameParameterScale;
	}

	public void setpMutateSameParameterScale(Double pMutateSameParameterScale) {
		this.pMutateSameParameterScale = pMutateSameParameterScale;
	}

	public String getPathToHRDScoresOutput() {
		return pathToHRDScoresOutput;
	}

	public void setPathToHRDScoresOutput(String pathToHRDScoresOutput) {
		this.pathToHRDScoresOutput = pathToHRDScoresOutput;
	}

	public Boolean doWriteHRDScoresToOutput() {
		return getPathToHRDScoresOutput() != null
				&& !getPathToHRDScoresOutput().equals("");
	}

	public String getPathToDomainWeightsDatabase() {
		return pathToDomainWeightsDatabase;
	}

	public void setPathToDomainWeightsDatabase(
			String pathToDomainWeightsDatabase) {
		this.pathToDomainWeightsDatabase = pathToDomainWeightsDatabase;
	}

	public String getPathToInterproResults4BlastHits() {
		return pathToInterproResults4BlastHits;
	}

	public void setPathToInterproResults4BlastHits(
			String pathToInterproResults4BlastHits) {
		this.pathToInterproResults4BlastHits = pathToInterproResults4BlastHits;
	}

	public String getPathToPfamResults4BlastHits() {
		return pathToPfamResults4BlastHits;
	}

	public void setPathToPfamResults4BlastHits(
			String pathToPfamResults4BlastHits) {
		this.pathToPfamResults4BlastHits = pathToPfamResults4BlastHits;
	}

	public Double getTokenScoreDomainSimilarityWeight() {
		return tokenScoreDomainSimilarityWeight;
	}

	public void setTokenScoreDomainSimilarityWeight(
			Double tokenScoreDomainSimilarityWeight) {
		this.tokenScoreDomainSimilarityWeight = tokenScoreDomainSimilarityWeight;
	}

	public Double getDescriptionScoreDomainSimilarityWeight() {
		return descriptionScoreDomainSimilarityWeight;
	}

	public void setDescriptionScoreDomainSimilarityWeight(
			Double descriptionScoreDomainSimilarityWeight) {
		this.descriptionScoreDomainSimilarityWeight = descriptionScoreDomainSimilarityWeight;
	}

	public String getComputeDomainSimilarityOn() {
		return computeDomainSimilarityOn;
	}

	public void setComputeDomainSimilarityOn(String computeDomainSimilarityOn) {
		if (computeDomainSimilarityOn != null)
			computeDomainSimilarityOn = computeDomainSimilarityOn.toLowerCase();
		this.computeDomainSimilarityOn = computeDomainSimilarityOn;
	}

	public boolean isDomainArchitectureSimilarityBasedOnPfamAnnotations() {
		return (getComputeDomainSimilarityOn() != null && getComputeDomainSimilarityOn()
				.equals("pfam"));
	}

	public boolean isWriteDomainArchitectureSimilarityScoresToOutput() {
		return writeDomainArchitectureSimilarityScoresToOutput;
	}

	public void setWriteDomainArchitectureSimilarityScoresToOutput(
			boolean writeDomainArchitectureSimilarityScoresToOutput) {
		this.writeDomainArchitectureSimilarityScoresToOutput = writeDomainArchitectureSimilarityScoresToOutput;
	}

	public Integer getDomainWeightTablePosition() {
		return domainWeightTablePosition;
	}

	public void setDomainWeightTablePosition(Integer domainWeightTablePosition) {
		this.domainWeightTablePosition = domainWeightTablePosition;
	}

}
