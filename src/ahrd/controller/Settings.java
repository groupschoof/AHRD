package ahrd.controller;

import static ahrd.controller.Utils.fromFile;
import static ahrd.controller.Utils.readFile;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Pattern;

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
	public static final String PROTEINS_FASTA_REGEX_KEY = "proteins_fasta_regex";
	public static final Pattern DEFAULT_PROTEINS_FASTA_REGEX = Pattern.compile("^>(?<accession>\\S+).*$");
	public static final String BLAST_DBS_KEY = "blast_dbs";
	public static final String BLAST_DB_WEIGHT_KEY = "weight";
	public static final String BLAST_RESULT_FILE_KEY = "file";
	public static final String BLAST_DATABASE_KEY = "database";
	public static final String BLAST_BLACKLIST_KEY = "blacklist";
	public static final String BLAST_FILTER_KEY = "filter";
	public static final String TOKEN_BLACKLIST_KEY = "token_blacklist";
	public static final String INTERPRO_DATABASE_KEY = "interpro_database";
	public static final String INTERPRO_RESULT_KEY = "interpro_result";
	public static final String GENE_ONTOLOGY_REFERENCE_KEY = "gene_ontology_reference";
	public static final String OUTPUT_KEY = "output";
	public static final String SIMULATED_ANNEALING_PATH_LOG_KEY = "path_log";
	public static final String WRITE_DESCIPTION_SUBSCORES_TO_OUTPUT = "write_description_subscores_to_output";
	public static final String WRITE_BEST_BLAST_HITS_TO_OUTPUT = "write_best_blast_hits_to_output";
	public static final String WRITE_TOKEN_SET_TO_OUTPUT = "write_token_set_to_output";
	public static final String HRD_SCORES_OUTPUT_PATH = "hrd_scores_output";
	public static final String TOKEN_SCORE_BIT_SCORE_WEIGHT = "token_score_bit_score_weight";
	public static final String TOKEN_SCORE_DATABASE_SCORE_WEIGHT = "token_score_database_score_weight";
	public static final String TOKEN_SCORE_OVERLAP_SCORE_WEIGHT = "token_score_overlap_score_weight";
	public static final String DESCRIPTION_SCORE_BIT_SCORE_WEIGHT = "description_score_bit_score_weight";
	public static final String GROUND_TRUTH_FASTA_KEY = "ground_truth_fasta";
	public static final String GROUND_TRUTH_DESCRIPTION_FILTER_KEY = "ground_truth_description_filter";
	public static final String GROUND_TRUTH_DESCRIPTION_BLACKLIST_KEY = "ground_truth_description_blacklist";
	public static final String GROUND_TRUTH_TOKEN_BLACKLIST_KEY = "ground_truth_token_blacklist";
	public static final String GROUND_TRUTH_FASTA_REGEX_KEY = "ground_truth_fasta_regex";
	public static final String F_MEASURE_BETA_PARAM_KEY = "f_measure_beta_parameter";
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
	public static final String SEQ_SIM_SEARCH_TABLE_COMMENT_LINE_REGEX_KEY = "seq_sim_search_table_comment_line_regex";
	public static final String SEQ_SIM_SEARCH_TABLE_SEP_KEY = "seq_sim_search_table_sep";
	public static final String SEQ_SIM_SEARCH_TABLE_QUERY_COL_KEY = "seq_sim_search_table_query_col";
	public static final String SEQ_SIM_SEARCH_TABLE_QUERY_COL_REGEX_KEY = "seq_sim_search_table_query_col_regex";
	public static final Pattern DEFAULT_SEQ_SIM_SEARCH_TABLE_QUERY_COL_REGEX = Pattern.compile("^(?<accession>.+)$");
	public static final String SEQ_SIM_SEARCH_TABLE_SUBJECT_COL_KEY = "seq_sim_search_table_subject_col";
	public static final String SEQ_SIM_SEARCH_TABLE_QUERY_START_COL_KEY = "seq_sim_search_table_query_start_col";
	public static final String SEQ_SIM_SEARCH_TABLE_QUERY_END_COL_KEY = "seq_sim_search_table_query_end_col";
	public static final String SEQ_SIM_SEARCH_TABLE_SUBJECT_START_COL_KEY = "seq_sim_search_table_subject_start_col";
	public static final String SEQ_SIM_SEARCH_TABLE_SUBJECT_END_COL_KEY = "seq_sim_search_table_subject_end_col";
	public static final String SEQ_SIM_SEARCH_TABLE_E_VALUE_COL_KEY = "seq_sim_search_table_e_value_col";
	public static final String SEQ_SIM_SEARCH_TABLE_BIT_SCORE_COL_KEY = "seq_sim_search_table_bit_score_col";
	public static final String FASTA_HEADER_REGEX_KEY = "fasta_header_regex";
	public static final Pattern DEFAULT_FASTA_HEADER_REGEX = Pattern
			.compile("^>(?<accession>\\S+)\\s+(?<description>.+?)(\\s+(((OS|os)=.+)|((GN|gn)=.+)))?$");
	public static final String SHORT_ACCESSION_REGEX_KEY = "short_accession_regex";
	public static final Pattern DEFAULT_SHORT_ACCESSION_REGEX = Pattern.compile("^[^|]+\\|(?<shortAccession>[^|]+)");
	public static final String GENE_ONTOLOGY_REFERENCE_REGEX_KEY = "gene_ontology_reference_regex";
	public static final Pattern DEFAULT_GENE_ONTOLOGY_REFERENCE_REGEX = Pattern
			.compile("^UniProtKB\\s+(?<shortAccession>\\S+)\\s+\\S+\\s+(?<goTerm>GO:\\d{7})");
	public static final String PREFER_REFERENCE_WITH_GO_ANNOS_KEY = "prefer_reference_with_go_annos";
	public static final String EVALUATE_ONLY_VALID_TOKENS_KEY = "evaluate_only_valid_tokens";
	public static final String DEFAULT_LINE_SEP = "(\r|\n)+";
	public static final String GO_DB_PATH_KEY = "go_db_path";
	public static final String GROUND_TRUTH_GO_ANNOTATIONS_PATH_KEY = "ground_truth_go_annotations";
	public static final String GO_F1_SIMPLE_KEY = "simple_GO_f1_scores";
	public static final String GO_F1_ANCESTRY_KEY = "ancestry_GO_f1_scores";
	public static final String GO_F1_SEMSIM_KEY = "semsim_GO_f1_scores";
	public static final String GO_SLIM_PATH_KEY = "go_slim";
	public static final Pattern GO_SLIM_FILE_GOTERM_REGEX = Pattern.compile("^id: (?<goTerm>GO:\\d{7})$");
	public static final String GENETIC_TRAINING_NUMBER_OF_GENERATIONS_KEY = "number_of_generations";
	public static final String GENETIC_TRAINING_POPULATION_SIZE_KEY = "population_size";
	public static final String COMPETITORS_KEY = "competitors";
	public static final String COMPETITOR_DESCRIPTIONS_FILE_KEY = "descriptions";
	public static final String COMPETITOR_GOA_FILE_KEY = "go_annotations";
	public static final String FIND_HIGHEST_POSSIBLE_GO_SCORE_KEY = "find_highest_possible_go_score";
	public static final String WRITE_FSCORE_DETAILS_TO_OUTPUT = "write_fscore_details_to_output";
	public static final String ACCESSION_GROUP_NAME = "accession";
	public static final String SHORT_ACCESSION_GROUP_NAME = "shortAccession";
	public static final String DESCRIPTION_GROUP_NAME = "description";
	public static final String GO_TERM_GROUP_NAME = "goTerm";
	public static final String FIND_HIGHEST_POSSIBLE_PRECISION_KEY = "find_highest_possible_precision";
	public static final String FIND_HIGHEST_POSSIBLE_RECALL_KEY = "find_highest_possible_recall";
	public static final String WRITE_EVALUATION_SUMMARY_KEY = "write_evaluation_summary";
	public static final String DESCRIPTION_SCORE_THRESHOLD = "description_score_threshold";
	
	/**
	 * Fields:
	 */
	private String pathToProteinsFasta;
	private Pattern proteinsFastaRegex;
	private String pathToGroundTruthFasta;
	private String pathToGroundTruthDescriptionBlacklist;
	private Set<String> groundTruthDescriptionBlacklist;
	private String pathToGroundTruthDescriptionFilter;
	private List<String> groundTruthDescriptionFilter;
	private String pathToGroundTruthTokenBlacklist;
	private Set<String> groundTruthTokenBlacklist = new HashSet<String>();
	private Pattern groundTruthFastaRegex;
	private String pathToInterproDatabase;
	private String pathToInterproResults;
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
	private Boolean writeBestBlastHitsToOutput = false;
	/**
	 * Forces AHRD to write out all internal scores (Sum(Token-Scores),
	 * Description- and Lexical-Scores, etc.
	 */
	private Boolean writeDescriptionSubScoresToOutputToOutput;
	/**
	 * F-Measure's Beta-Parameter as set in the input.yml or default 1.0
	 */
	private Double fMeasureBetaParameter = 1.0;
	private Map<String, Map<String, String>> blastDbSettings = new HashMap<String, Map<String, String>>();
	private List<String> sortedBlastDatabaseNames;
	private Map<String, Set<String>> blastResultsBlacklists = new HashMap<String, Set<String>>();
	private Map<String, List<String>> blastResultsFilter = new HashMap<String, List<String>>();
	private Map<String, Set<String>> tokenBlacklists = new HashMap<String, Set<String>>();
	/**
	 * For the <strong>simulated annealing</strong> algorithm, this will be
	 * current temperature. (Default is 75000)
	 */
	private Integer temperature = 75000;
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
	private Double optimizationAcceptanceProbabilityScalingFactor = 2500000000.0;
	/**
	 * In simulated annealing optimization each cycle has to mutate the current
	 * Parameter-Set to generate a neighboring set in parameter space. The
	 * random value used to add or subtract to a single parameter is Gaussian
	 * distributed and has the following mean:
	 */
	private Double mutatorMean = 0.25;
	/**
	 * In simulated annealing optimization each cycle has to mutate the current
	 * Parameter-Set to generate a neighboring set in parameter space. The
	 * random value used to add or subtract to a single parameter is Gaussian
	 * distributed and has the following standard deviation:
	 */
	private Double mutatorDeviation = 0.15;
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
	 * The following fields control how the result table of a sequence
	 * similarity search is parsed. All concerned fields start with
	 * 'seqSimSearchTable'.
	 */
	private Pattern seqSimSearchTableCommentLineRegex = null;
	private String seqSimSearchTableSep = "\t";
	private Integer seqSimSearchTableQueryCol = 0;
	private Pattern seqSimSearchTableQueryColRegex;
	private Integer seqSimSearchTableSubjectCol = 1;
	private Integer seqSimSearchTableQueryStartCol = 6;
	private Integer seqSimSearchTableQueryEndCol = 7;
	private Integer seqSimSearchTableSubjectStartCol = 8;
	private Integer seqSimSearchTableSubjectEndCol = 9;
	private Integer seqSimSearchTableEValueCol = 10;
	private Integer seqSimSearchTableBitScoreCol = 11;
	/**
	 * If set to true AHRD will choose the highest scoring BlastResult WITH GO
	 * Annotations as donor for a query protein's HRD. If no BlastResult has GO
	 * annotations AHRD works "as normal".
	 */
	private Boolean preferReferenceWithGoAnnos = true;
	/**
	 * If set to TRUE the AHRD Evaluation Score is based ONLY on tokens that
	 * pass the Blacklisting. Otherwise all Tokens are submitted to evaluation.
	 */
	private Boolean evaluateOnlyValidTokens = true;
	/**
	 * The path in witch to keep: - Downloaded reviewed part of Uniprot -
	 * Downloaded Gene Ontology mysql database dump - Serialized copy of the
	 * GOdatabase, generated when it was needed the for first time
	 */
	private String pathToGoDatabase;
	/**
	 * Path to file containing GO annotations of the query proteins.
	 * One of the requirements to trigger the evaluation of AHRDs GO annotations
	 */
	private String pathToGroundTruthGoAnnotations;
	/**
	 * Toggle set to induce the calculation and output of an F1-score based on
	 * ground truth and prediction GO annotations alone. Is set to true as default
	 * GO F1-score if non of the three GO F1-score flags are toggled.
	 */
	private Boolean calculateSimpleGoF1Scores = false;
	/**
	 * Toggle set to induce the calculation and output of an F1-score based on
	 * ground truth and prediction GO annotations extended to their complete
	 * ancestry.
	 */
	private Boolean calculateAncestryGoF1Scores = false;
	/**
	 * Toggle set to induce the calculation and output of an F1-score based on
	 * the semantic similarity (based on term information content) of ground truth
	 * and prediction GO annotations.
	 */
	private Boolean calculateSemSimGoF1Scores = false;
	/**
	 * Path to file containing a subset of the Gene Ontology. Triggers the
	 * annotation of the query proteins with this broad categorical GO terms
	 * instead of the detailed ones determined from the blast results.
	 * Useful to give a summary of the GO annotation of for example a genome, microarray or cDNA collection.
	 */
	private String pathToGoSlimFile;
	/**
	 * The number of generation to be consecutively evolved and evaluated when performing parameter optimization in the genetic trainer
	 */
	private int numberOfGenerations = 100;
	/**
	 * The size of each generation to be evolved and evaluated when performing parameter optimization in the genetic trainer
	 */
	private int populationSize = 200;
	/**
	 * Competitors to be compared to AHRD in evaluation run
	 */
	private Map<String, Map<String, String>> competitorSettings = new HashMap<String, Map<String, String>>();
	/**
	 * Evaluation or Optimization might be interested in the highest possibly
	 * achievable score for go annotations:
	 */
	private boolean findHighestPossibleGoScore = false;
	/**
	 * Adds the precision and recall to the output of all (go based) F-scores 
	 */
	private boolean writeFscoreDetailsToOutput = false;

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
	 * Default description line blacklist.
	 * Is applied to all blast databases that don't have a description line blacklist specified.
	 */
	private Set<String> defaultBlastResultsBlacklist = new HashSet<String>();
	/**
	 * Default description line filter.
	 * Is applied to all blast databases that don't have a description line filter specified.
	 */
	private List<String> defaultBlastResultsFilter = new ArrayList<String>();
	/**
	 * Default token blacklist.
	 * Is applied to all blast databases that don't have a token blacklist specified.
	 */
	private Set<String> defaultTokenBlacklist = new HashSet<String>();
	/**
	 * Triggers the calculation and output of the highest possible precision of descriptions and gene ontology terms.
	 */
	private boolean findHighestPossiblePrecision = false;
	/**
	 * Triggers the calculation and output of the highest possible recall of descriptions and gene ontology terms.
	 */
	private boolean findHighestPossibleRecall = false;
	/**
	 * Triggers the output of averages and coverages of all scores at the end of the evalution output
	 */
	private boolean writeEvaluationSummary = false;

	/**
	 * Initializes an Instance with content read from a YML-File:
	 * 
	 * @throws IOException
	 */
	@SuppressWarnings("unchecked")
	public void initialize(String pathToYml) throws IOException {
		YamlReader reader = new YamlReader(new FileReader(pathToYml));
		Map<String, Object> input = (Map<String, Object>) reader.read();
		this.blastDbSettings = (Map<String, Map<String, String>>) input.get(BLAST_DBS_KEY);
		this.setPathToProteinsFasta((String) input.get(PROTEINS_FASTA_KEY));
		if (input.get(PROTEINS_FASTA_REGEX_KEY) != null) {
			this.setProteinsFastaRegex(Pattern.compile((String) input.get(PROTEINS_FASTA_REGEX_KEY)));
		} else {
			this.setProteinsFastaRegex(DEFAULT_PROTEINS_FASTA_REGEX);
		}
		// If started to train the algorithm ground truth descriptions are stored in this file:
		setPathToGroundTruthFasta((String) input.get(GROUND_TRUTH_FASTA_KEY));
		this.setPathToInterproDatabase((String) input.get(INTERPRO_DATABASE_KEY));
		this.setPathToInterproResults((String) input.get(INTERPRO_RESULT_KEY));
		this.setPathToOutput((String) input.get(OUTPUT_KEY));
		if (input.get(HRD_SCORES_OUTPUT_PATH) != null && !input.get(HRD_SCORES_OUTPUT_PATH).equals(""))
			setPathToHRDScoresOutput((String) input.get(HRD_SCORES_OUTPUT_PATH));
		// Trainer logs path through parameter-space here:
		if (input.get(SIMULATED_ANNEALING_PATH_LOG_KEY) != null)
			setPathToSimulatedAnnealingPathLog((String) input.get(SIMULATED_ANNEALING_PATH_LOG_KEY));
		this.setTokenScoreBitScoreWeight(Double.parseDouble((String) input.get(TOKEN_SCORE_BIT_SCORE_WEIGHT)));
		this.setTokenScoreDatabaseScoreWeight(Double.parseDouble((String) input.get(TOKEN_SCORE_DATABASE_SCORE_WEIGHT)));
		this.setTokenScoreOverlapScoreWeight(Double.parseDouble((String) input.get(TOKEN_SCORE_OVERLAP_SCORE_WEIGHT)));
		this.setWriteTokenSetToOutput(Boolean.parseBoolean((String) input.get(WRITE_TOKEN_SET_TO_OUTPUT)));
		// Writing best blast hits to output is only supported in evaluation mode. So otherwise the default 'false' will be kept.
		if (isInEvaluationMode()) {
			this.setWriteBestBlastHitsToOutput(Boolean.parseBoolean((String) input.get(WRITE_BEST_BLAST_HITS_TO_OUTPUT)));
		}
		this.setWriteDesciptionSubScoresToOutput(Boolean.parseBoolean((String) input.get(WRITE_DESCIPTION_SUBSCORES_TO_OUTPUT)));
		this.setOutputFasta(Boolean.parseBoolean((String) input.get(OUTPUT_FASTA_KEY)));
		if (input.get(BLAST_BLACKLIST_KEY) != null) {
			this.setDefaultBlastResultsBlacklist(new HashSet<String>(fromFile((String) input.get(BLAST_BLACKLIST_KEY))));
		}
		if (input.get(BLAST_FILTER_KEY) != null) {
			this.setDefaultBlastResultsFilter(fromFile((String) input.get(BLAST_FILTER_KEY)));
		}
		if (input.get(TOKEN_BLACKLIST_KEY) != null) {
			this.setDefaultTokenBlacklist(new HashSet<String>(fromFile((String) input.get(TOKEN_BLACKLIST_KEY))));
		}
		// Generate the Blacklists and Filters for each Blast-Database from their appropriate files:
		for (String blastDatabaseName : getBlastDatabases()) {
			if (getPathToBlastResultsBlacklist(blastDatabaseName) != null) {
				this.blastResultsBlacklists.put(blastDatabaseName, new HashSet<String>(fromFile(getPathToBlastResultsBlacklist(blastDatabaseName))));
			}
			if (getPathToBlastResultsFilter(blastDatabaseName) != null) {
				this.blastResultsFilter.put(blastDatabaseName, fromFile(getPathToBlastResultsFilter(blastDatabaseName)));
			}
			if (getPathToTokenBlacklist(blastDatabaseName) != null) {
				this.tokenBlacklists.put(blastDatabaseName, new HashSet<String>(fromFile(getPathToTokenBlacklist(blastDatabaseName))));
			}
			// Set Database-Weights and Description-Score-Bit-Score-Weight:
			this.getParameters().setBlastDbWeight(blastDatabaseName,
					this.getBlastDbSettings(blastDatabaseName).get(Settings.BLAST_DB_WEIGHT_KEY));
			this.getParameters().setDescriptionScoreBitScoreWeight(blastDatabaseName,
					this.getBlastDbSettings(blastDatabaseName).get(Settings.DESCRIPTION_SCORE_BIT_SCORE_WEIGHT));
		}
		// If started in evaluation-mode the F-Measure's Beta-Parameter can be set
		// to some other value than 1.0
		if (input.get(F_MEASURE_BETA_PARAM_KEY) != null)
			this.fMeasureBetaParameter = Double.parseDouble((String) input.get(F_MEASURE_BETA_PARAM_KEY));
		// Simulated Annealing can be started with custom temperature and value
		// it is cooled-down by each step:
		if (input.get(TEMPERATURE_KEY) != null)
			this.setTemperature(Integer.parseInt((String) input.get(TEMPERATURE_KEY)));
		if (input.get(COOL_DOWN_BY_KEY) != null)
			this.coolDownBy = Integer.parseInt((String) input.get(COOL_DOWN_BY_KEY));
		if (input.get(OPTIMIZATION_ACCEPTANCE_PROBABILITY_SCALING_FACTOR_KEY) != null)
			this.setOptimizationAcceptanceProbabilityScalingFactor(
					Double.parseDouble((String) input.get(OPTIMIZATION_ACCEPTANCE_PROBABILITY_SCALING_FACTOR_KEY)));
		if (input.get(MUTATOR_MEAN_KEY) != null)
			this.setMutatorMean(Double.parseDouble((String) input.get(MUTATOR_MEAN_KEY)));
		if (input.get(MUTATOR_DEVIATION_KEY) != null)
			this.setMutatorDeviation(Double.parseDouble((String) input.get(MUTATOR_DEVIATION_KEY)));
		this.setRememberSimulatedAnnealingPath(Boolean.parseBoolean((String) input.get(REMEMBER_SIMULATED_ANNEALING_PATH_KEY)));
		if (input.get(P_MUTATE_SAME_PARAMETER_SCALE_KEY) != null)
			this.setpMutateSameParameterScale(Double.parseDouble((String) input.get(P_MUTATE_SAME_PARAMETER_SCALE_KEY)));
		// Evaluation or Optimization might be interested in the highest
		// possibly achievable evaluation-score:
		this.setFindHighestPossibleEvaluationScore(Boolean.parseBoolean((String) input.get(FIND_HIGHEST_POSSIBLE_EVALUATION_SCORE_KEY)));
		// Set any non default parameters controlling, how sequence similarity
		// search result tables are parsed:
		if (input.get(SEQ_SIM_SEARCH_TABLE_COMMENT_LINE_REGEX_KEY) != null) {
			this.setSeqSimSearchTableCommentLineRegex(
					Pattern.compile(input.get(SEQ_SIM_SEARCH_TABLE_COMMENT_LINE_REGEX_KEY).toString()));
		}
		if (input.get(SEQ_SIM_SEARCH_TABLE_SEP_KEY) != null) {
			this.setSeqSimSearchTableSep(input.get(SEQ_SIM_SEARCH_TABLE_SEP_KEY).toString());
		}
		if (input.get(SEQ_SIM_SEARCH_TABLE_QUERY_COL_KEY) != null) {
			this.setSeqSimSearchTableQueryCol(Integer.parseInt(input.get(SEQ_SIM_SEARCH_TABLE_QUERY_COL_KEY).toString()));
		}
		if (input.get(SEQ_SIM_SEARCH_TABLE_QUERY_COL_REGEX_KEY) != null) {
			this.setSeqSimSearchTableQueryColRegex(Pattern.compile((String) input.get(SEQ_SIM_SEARCH_TABLE_QUERY_COL_REGEX_KEY)));
		} else {
			this.setSeqSimSearchTableQueryColRegex(DEFAULT_SEQ_SIM_SEARCH_TABLE_QUERY_COL_REGEX);
		}
		if (input.get(SEQ_SIM_SEARCH_TABLE_SUBJECT_COL_KEY) != null) {
			this.setSeqSimSearchTableSubjectCol(
					Integer.parseInt(input.get(SEQ_SIM_SEARCH_TABLE_SUBJECT_COL_KEY).toString()));
		}
		if (input.get(SEQ_SIM_SEARCH_TABLE_QUERY_START_COL_KEY) != null) {
			this.setSeqSimSearchTableQueryStartCol(
					Integer.parseInt(input.get(SEQ_SIM_SEARCH_TABLE_QUERY_START_COL_KEY).toString()));
		}
		if (input.get(SEQ_SIM_SEARCH_TABLE_QUERY_END_COL_KEY) != null) {
			this.setSeqSimSearchTableQueryEndCol(
					Integer.parseInt(input.get(SEQ_SIM_SEARCH_TABLE_QUERY_END_COL_KEY).toString()));
		}
		if (input.get(SEQ_SIM_SEARCH_TABLE_SUBJECT_START_COL_KEY) != null) {
			this.setSeqSimSearchTableSubjectStartCol(
					Integer.parseInt(input.get(SEQ_SIM_SEARCH_TABLE_SUBJECT_START_COL_KEY).toString()));
		}
		if (input.get(SEQ_SIM_SEARCH_TABLE_SUBJECT_END_COL_KEY) != null) {
			this.setSeqSimSearchTableSubjectEndCol(
					Integer.parseInt(input.get(SEQ_SIM_SEARCH_TABLE_SUBJECT_END_COL_KEY).toString()));
		}
		if (input.get(SEQ_SIM_SEARCH_TABLE_E_VALUE_COL_KEY) != null) {
			this.setSeqSimSearchTableEValueCol(Integer.parseInt(input.get(SEQ_SIM_SEARCH_TABLE_E_VALUE_COL_KEY).toString()));
		}
		if (input.get(SEQ_SIM_SEARCH_TABLE_BIT_SCORE_COL_KEY) != null) {
			this.setSeqSimSearchTableBitScoreCol(
					Integer.parseInt(input.get(SEQ_SIM_SEARCH_TABLE_BIT_SCORE_COL_KEY).toString()));
		}
		if (input.get(PREFER_REFERENCE_WITH_GO_ANNOS_KEY) != null) { 
			this.setPreferReferenceWithGoAnnos(Boolean.parseBoolean((String) input.get(PREFER_REFERENCE_WITH_GO_ANNOS_KEY)));
		}
		if (input.get(EVALUATE_ONLY_VALID_TOKENS_KEY) != null) {
			this.setEvaluateOnlyValidTokens(Boolean.parseBoolean((String) input.get(EVALUATE_ONLY_VALID_TOKENS_KEY)));
		}
		if (input.get(GROUND_TRUTH_DESCRIPTION_BLACKLIST_KEY) != null) {
			this.setPathToGroundTruthDescriptionBlacklist(input.get(GROUND_TRUTH_DESCRIPTION_BLACKLIST_KEY).toString());
			this.setGroundTruthDescriptionBlacklist(new HashSet<String>(fromFile(getPathToGroundTruthDescriptionBlacklist())));
		}
		if (input.get(GROUND_TRUTH_DESCRIPTION_FILTER_KEY) != null) {
			this.setPathToGroundTruthDescriptionFilter(input.get(GROUND_TRUTH_DESCRIPTION_FILTER_KEY).toString());
			this.setGroundTruthDescriptionFilter(fromFile(getPathToGroundTruthDescriptionFilter()));
		}
		// If the ground_truth_token_blacklist parameter is set the provided file will be used as blacklist for the ground truth tokens.
		// Otherwise the evaluate_only_valid_tokens parameter is checked. If set to 'true' (its default) the defaultTokenBlacklist is used to filter the ground truth tokens.
		// If evaluate_only_valid_tokens is set to 'false' the ground truth token black list remains empty,
		// which will result in ALL ground truth tokens being used in the evaluation (so no filtering takes place).
		if (input.get(GROUND_TRUTH_TOKEN_BLACKLIST_KEY) != null) {
			this.setPathToGroundTruthTokenBlacklist(input.get(GROUND_TRUTH_TOKEN_BLACKLIST_KEY).toString());
			this.setGroundTruthTokenBlacklist(new HashSet<String>(fromFile(getPathToGroundTruthTokenBlacklist())));
		} else {
			if (this.getEvaluateOnlyValidTokens()) {
				this.setGroundTruthTokenBlacklist(this.getDefaultTokenBlacklist());
			}
		}
		if (input.get(GROUND_TRUTH_FASTA_REGEX_KEY) != null) {
			this.setGroundTruthFastaRegex(Pattern.compile((String) input.get(GROUND_TRUTH_FASTA_REGEX_KEY)));
		} else {
			this.setGroundTruthFastaRegex(DEFAULT_FASTA_HEADER_REGEX);
		}
		if (input.get(GO_DB_PATH_KEY) != null) {
			this.setPathToGoDatabase(input.get(GO_DB_PATH_KEY).toString());
		}
		if (input.get(GROUND_TRUTH_GO_ANNOTATIONS_PATH_KEY) != null) {
			this.setPathToGroundTruthGoAnnotations(input.get(GROUND_TRUTH_GO_ANNOTATIONS_PATH_KEY).toString());
		}
		this.setCalculateAncestryGoF1Scores(Boolean.parseBoolean((String) input.get(GO_F1_ANCESTRY_KEY)));
		this.setCalculateSemSimGoF1Scores(Boolean.parseBoolean((String) input.get(GO_F1_SEMSIM_KEY)));
		// If non of the GO F1 Keys is specified in the YML input the simple
		// version is used as default
		if (input.get(GO_F1_SIMPLE_KEY) != null) {
			this.setCalculateSimpleGoF1Scores(Boolean.parseBoolean((String) input.get(GO_F1_SIMPLE_KEY)));
		} else {
			if (input.get(GO_F1_ANCESTRY_KEY) == null && input.get(GO_F1_SEMSIM_KEY) == null) {
				this.setCalculateSimpleGoF1Scores(true);
			}
		}
		if (input.get(GO_SLIM_PATH_KEY) != null) {
			this.setPathToGoSlimFile(input.get(GO_SLIM_PATH_KEY).toString());
		}
		if (input.get(GENETIC_TRAINING_NUMBER_OF_GENERATIONS_KEY) != null) {
			this.setNumberOfGenerations(Integer.parseInt((String) input.get(GENETIC_TRAINING_NUMBER_OF_GENERATIONS_KEY)));
		}
		if (input.get(GENETIC_TRAINING_POPULATION_SIZE_KEY) != null) {
			this.setPopulationSize(Integer.parseInt((String) input.get(GENETIC_TRAINING_POPULATION_SIZE_KEY)));
		}
		if (input.get(COMPETITORS_KEY) != null) {
			setCompetitorSettings((Map<String, Map<String, String>>) input.get(COMPETITORS_KEY));
		}
		this.setFindHighestPossibleGoScore(Boolean.parseBoolean((String) input.get(FIND_HIGHEST_POSSIBLE_GO_SCORE_KEY)));
		this.setWriteFscoreDetailsToOutput(Boolean.parseBoolean((String) input.get(WRITE_FSCORE_DETAILS_TO_OUTPUT)));
		this.setFindHighestPossiblePrecision(Boolean.parseBoolean((String) input.get(FIND_HIGHEST_POSSIBLE_PRECISION_KEY)));
		this.setFindHighestPossibleRecall(Boolean.parseBoolean((String) input.get(FIND_HIGHEST_POSSIBLE_RECALL_KEY)));
		this.setWriteEvaluationSummary(Boolean.parseBoolean((String) input.get(WRITE_EVALUATION_SUMMARY_KEY)));
		if (input.get(DESCRIPTION_SCORE_THRESHOLD) != null) {
			this.setDescriptionScoreThreshold(Double.parseDouble((String) input.get(DESCRIPTION_SCORE_THRESHOLD)));
		}
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
		if (getPathToInterproDatabase() == null || getPathToInterproResults() == null)
			return false;
		// ELSE:
		File iprDb = new File(getPathToInterproDatabase());
		File iprRes = new File(getPathToInterproResults());
		return (iprDb.canRead() && iprDb.length() > 0 && iprRes.canRead() && iprRes.length() > 0);
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
	
	public void setRememberSimulatedAnnealingPath(boolean rememberSimulatedAnnealingPath) {
		this.rememberSimulatedAnnealingPath = rememberSimulatedAnnealingPath;
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
			this.sortedBlastDatabaseNames = new ArrayList<String>(getBlastDatabases());
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
		return getParameters().getDescriptionScoreBitScoreWeight(blastDatabaseName);
	}

	public void setDescriptionScoreBitScoreWeight(String blastDatabaseName, String dsbsw) {
		getParameters().setDescriptionScoreBitScoreWeight(blastDatabaseName, dsbsw);
	}

	public String getPathToBlastResults(String blastDatabaseName) {
		return getBlastDbSettings(blastDatabaseName).get(BLAST_RESULT_FILE_KEY);
	}

	public String getPathToBlastDatabase(String blastDatabaseName) {
		return getBlastDbSettings(blastDatabaseName).get(BLAST_DATABASE_KEY);
	}

	public Pattern getFastaHeaderRegex(String blastDatabaseName) {
		return (getBlastDbSettings(blastDatabaseName).containsKey(FASTA_HEADER_REGEX_KEY))
				? Pattern.compile(getBlastDbSettings(blastDatabaseName).get(FASTA_HEADER_REGEX_KEY).toString())
				: DEFAULT_FASTA_HEADER_REGEX;
	}

	public Pattern getShortAccessionRegex(String blastDatabaseName) {
		return (getBlastDbSettings(blastDatabaseName).containsKey(SHORT_ACCESSION_REGEX_KEY))
				? Pattern.compile(getBlastDbSettings(blastDatabaseName).get(SHORT_ACCESSION_REGEX_KEY).toString())
				: DEFAULT_SHORT_ACCESSION_REGEX;
	}

	private String getPathToBlastResultsBlacklist(String blastDatabaseName) {
		return getBlastDbSettings(blastDatabaseName).get(BLAST_BLACKLIST_KEY);
	}

	public Set<String> getBlastResultsBlacklist() {
		return this.getDefaultBlastResultsBlacklist();
	}
	
	public Set<String> getBlastResultsBlacklist(String blastDatabaseName) {
		if (this.blastResultsBlacklists.get(blastDatabaseName) != null) {
			return this.blastResultsBlacklists.get(blastDatabaseName);
		}
		return this.getDefaultBlastResultsBlacklist();
	}

	private String getPathToBlastResultsFilter(String blastDatabaseName) {
		return getBlastDbSettings(blastDatabaseName).get(BLAST_FILTER_KEY);
	}

	public List<String> getBlastResultsFilter() {
		return this.getDefaultBlastResultsFilter();
	}
	
	public List<String> getBlastResultsFilter(String blastDatabaseName) {
		if (this.blastResultsFilter.get(blastDatabaseName) != null) {
			return this.blastResultsFilter.get(blastDatabaseName);
		}
		return this.getDefaultBlastResultsFilter();
	}

	private String getPathToTokenBlacklist(String blastDatabaseName) {
		return getBlastDbSettings(blastDatabaseName).get(TOKEN_BLACKLIST_KEY);
	}

	public Set<String> getTokenBlacklist() {
		return this.getDefaultTokenBlacklist();
	}
	
	public Set<String> getTokenBlacklist(String blastDatabaseName) {
		if (this.tokenBlacklists.get(blastDatabaseName) != null) {
			return this.tokenBlacklists.get(blastDatabaseName);	
		}
		return this.getDefaultTokenBlacklist();
	}

	private String getPathToProteinsFasta() {
		return pathToProteinsFasta;
	}

	public String getProteinsFasta() throws IOException {
		return readFile(getPathToProteinsFasta());
	}

	private String getPathToGroundTruthFasta() {
		return pathToGroundTruthFasta;
	}

	public String getGroundTruthFasta() throws IOException {
		return readFile(getPathToGroundTruthFasta());
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
		return getPathToInterproDatabase() != null && (new File(getPathToInterproDatabase())).exists()
				&& getPathToInterproResults() != null && (new File(getPathToInterproResults())).exists();
	}

	public void setPathToInterproResults(String pathToInterproResults) {
		this.pathToInterproResults = pathToInterproResults;
	}

	public String getPathToGeneOntologyReference(String blastDatabaseName) {
		return getBlastDbSettings(blastDatabaseName).get(GENE_ONTOLOGY_REFERENCE_KEY);
	}

	public boolean hasGeneOntologyAnnotations() {
		boolean result = false;
		for (String blastDatabaseName : getBlastDatabases()) {
			if (getPathToGeneOntologyReference(blastDatabaseName) != null
					&& new File(getPathToGeneOntologyReference(blastDatabaseName)).exists()) {
				result = true;
			}
		}
		return result;
	}

	public boolean hasGeneOntologyAnnotation(String blastDatabaseName) {
		if (getPathToGeneOntologyReference(blastDatabaseName) != null
				&& new File(getPathToGeneOntologyReference(blastDatabaseName)).exists()) {
			return true;
		}
		return false;
	}

	public void setPathToGeneOntologyReference(String blastDatabaseName, String path) {
		getBlastDbSettings(blastDatabaseName).put(GENE_ONTOLOGY_REFERENCE_KEY, path);
	}

	public void removeAllPathToGeneOntologyReferemces() {
		for (String blastDatabaseName : getBlastDatabases()) {
			getBlastDbSettings(blastDatabaseName).remove(GENE_ONTOLOGY_REFERENCE_KEY);
		}
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
		this.getParameters().setTokenScoreBitScoreWeight(tokenScoreBitScoreWeight);
	}

	public Double getTokenScoreDatabaseScoreWeight() {
		return getParameters().getTokenScoreDatabaseScoreWeight();
	}

	public void setTokenScoreDatabaseScoreWeight(Double tokenScoreDatabaseScoreWeight) {
		this.getParameters().setTokenScoreDatabaseScoreWeight(tokenScoreDatabaseScoreWeight);
	}

	public Double getTokenScoreOverlapScoreWeight() {
		return getParameters().getTokenScoreOverlapScoreWeight();
	}

	public void setTokenScoreOverlapScoreWeight(Double tokenScoreOverlapScoreWeight) {
		this.getParameters().setTokenScoreOverlapScoreWeight(tokenScoreOverlapScoreWeight);
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

	public Boolean getWriteDescriptionSubScoresToOutput() {
		return writeDescriptionSubScoresToOutputToOutput;
	}

	public void setWriteDesciptionSubScoresToOutput(Boolean writeScoresToOutput) {
		this.writeDescriptionSubScoresToOutputToOutput = writeScoresToOutput;
	}

	public Map<String, Map<String, String>> getBlastDbSettings() {
		return blastDbSettings;
	}

	public void setBlastDbSettings(Map<String, Map<String, String>> blastDbSettings) {
		this.blastDbSettings = blastDbSettings;
	}

	public void setPathToGroundTruthFasta(String pathToGroundTruthFasta) {
		this.pathToGroundTruthFasta = pathToGroundTruthFasta;
	}

	public boolean isInEvaluationMode() {
		return (getPathToGroundTruthFasta() != null && getPathToGroundTruthFasta() != "");
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


	public Double getAvgTrainingScore() {
		return getParameters().getAvgTrainingScore();
	}

	public void setAvgTrainingScore(Double avgEvaluationScore) {
		this.getParameters().setAvgTrainingScore(avgEvaluationScore);
	}
	
	public Double getAvgPrecision() {
		return getParameters().getAvgPrecision();
	}

	public void setAvgPrecision(Double avgPrecision) {
		this.getParameters().setAvgPrecision(avgPrecision);
	}

	public Double getAvgRecall() {
		return getParameters().getAvgRecall();
	}

	public void setAvgRecall(Double avgRecall) {
		this.getParameters().setAvgRecall(avgRecall);
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

	public void setPathToSimulatedAnnealingPathLog(String pathToSimulatedAnnealingPathLog) {
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
		return getPathToHRDScoresOutput() != null && !getPathToHRDScoresOutput().equals("");
	}

	public Integer getSeqSimSearchTableQueryCol() {
		return seqSimSearchTableQueryCol;
	}

	public void setSeqSimSearchTableQueryCol(Integer seqSimSearchTableQueryCol) {
		this.seqSimSearchTableQueryCol = seqSimSearchTableQueryCol;
	}

	public Integer getSeqSimSearchTableSubjectCol() {
		return seqSimSearchTableSubjectCol;
	}

	public void setSeqSimSearchTableSubjectCol(Integer seqSimSearchTableSubjectCol) {
		this.seqSimSearchTableSubjectCol = seqSimSearchTableSubjectCol;
	}

	public Integer getSeqSimSearchTableQueryStartCol() {
		return seqSimSearchTableQueryStartCol;
	}

	public void setSeqSimSearchTableQueryStartCol(Integer seqSimSearchTableQueryStartCol) {
		this.seqSimSearchTableQueryStartCol = seqSimSearchTableQueryStartCol;
	}

	public Integer getSeqSimSearchTableQueryEndCol() {
		return seqSimSearchTableQueryEndCol;
	}

	public void setSeqSimSearchTableQueryEndCol(Integer seqSimSearchTableQueryEndCol) {
		this.seqSimSearchTableQueryEndCol = seqSimSearchTableQueryEndCol;
	}

	public Integer getSeqSimSearchTableSubjectStartCol() {
		return seqSimSearchTableSubjectStartCol;
	}

	public void setSeqSimSearchTableSubjectStartCol(Integer seqSimSearchTableSubjectStartCol) {
		this.seqSimSearchTableSubjectStartCol = seqSimSearchTableSubjectStartCol;
	}

	public Integer getSeqSimSearchTableSubjectEndCol() {
		return seqSimSearchTableSubjectEndCol;
	}

	public void setSeqSimSearchTableSubjectEndCol(Integer seqSimSearchTableSubjectEndCol) {
		this.seqSimSearchTableSubjectEndCol = seqSimSearchTableSubjectEndCol;
	}

	public Pattern getSeqSimSearchTableCommentLineRegex() {
		return seqSimSearchTableCommentLineRegex;
	}

	public void setSeqSimSearchTableCommentLineRegex(Pattern seqSimSearchTableCommentLineRegex) {
		this.seqSimSearchTableCommentLineRegex = seqSimSearchTableCommentLineRegex;
	}

	public String getSeqSimSearchTableSep() {
		return seqSimSearchTableSep;
	}

	public void setSeqSimSearchTableSep(String seqSimSearchTableSep) {
		this.seqSimSearchTableSep = seqSimSearchTableSep;
	}

	public Integer getSeqSimSearchTableEValueCol() {
		return seqSimSearchTableEValueCol;
	}

	public void setSeqSimSearchTableEValueCol(Integer seqSimSearchTableEValueCol) {
		this.seqSimSearchTableEValueCol = seqSimSearchTableEValueCol;
	}

	public Integer getSeqSimSearchTableBitScoreCol() {
		return seqSimSearchTableBitScoreCol;
	}

	public void setSeqSimSearchTableBitScoreCol(Integer seqSimSearchTableBitScoreCol) {
		this.seqSimSearchTableBitScoreCol = seqSimSearchTableBitScoreCol;
	}

	/**
	 * Either returns the custom regular expression pattern used to parse the
	 * provided Gene Ontology annotion (GOA) reference or returns the default
	 * pattern designed to work for UniprotKB GOA files.
	 * 
	 * @return Pattern
	 */
	public Pattern getGoReferenceRegex(String blastDatabaseName) {
		if (getBlastDbSettings(blastDatabaseName).get(GENE_ONTOLOGY_REFERENCE_REGEX_KEY) != null) {
			return Pattern.compile(getBlastDbSettings(blastDatabaseName).get(GENE_ONTOLOGY_REFERENCE_REGEX_KEY).toString());
		}
		return DEFAULT_GENE_ONTOLOGY_REFERENCE_REGEX;
	}

	public Boolean getPreferReferenceWithGoAnnos() {
		return preferReferenceWithGoAnnos;
	}

	public void setPreferReferenceWithGoAnnos(Boolean preferReferenceWithGoAnnos) {
		this.preferReferenceWithGoAnnos = preferReferenceWithGoAnnos;
	}

	public Boolean getEvaluateOnlyValidTokens() {
		return evaluateOnlyValidTokens;
	}

	public void setEvaluateOnlyValidTokens(Boolean evaluateValidTokens) {
		this.evaluateOnlyValidTokens = evaluateValidTokens;
	}

	public String getPathToGroundTruthDescriptionFilter() {
		return pathToGroundTruthDescriptionFilter;
	}

	public void setPathToGroundTruthDescriptionFilter(String pathToGroundTruthDescriptionFilter) {
		this.pathToGroundTruthDescriptionFilter = pathToGroundTruthDescriptionFilter;
	}

	public String getPathToGroundTruthDescriptionBlacklist() {
		return pathToGroundTruthDescriptionBlacklist;
	}

	public void setPathToGroundTruthDescriptionBlacklist(String pathToGroundTruthDescriptionBlacklist) {
		this.pathToGroundTruthDescriptionBlacklist = pathToGroundTruthDescriptionBlacklist;
	}

	public String getPathToGroundTruthTokenBlacklist() {
		return pathToGroundTruthTokenBlacklist;
	}

	public void setPathToGroundTruthTokenBlacklist(String pathToGroundTruthTokenBlacklist) {
		this.pathToGroundTruthTokenBlacklist = pathToGroundTruthTokenBlacklist;
	}

	public Set<String> getGroundTruthDescriptionBlacklist() {
		return groundTruthDescriptionBlacklist;
	}

	public void setGroundTruthDescriptionBlacklist(Set<String> groundTruthDescriptionBlacklist) {
		this.groundTruthDescriptionBlacklist = groundTruthDescriptionBlacklist;
	}

	public List<String> getGroundTruthDescriptionFilter() {
		return groundTruthDescriptionFilter;
	}

	public void setGroundTruthDescriptionFilter(List<String> groundTruthDescriptionFilter) {
		this.groundTruthDescriptionFilter = groundTruthDescriptionFilter;
	}

	public Set<String> getGroundTruthTokenBlacklist() {
		return groundTruthTokenBlacklist;
	}

	public void setGroundTruthTokenBlacklist(Set<String> groundTruthTokenBlacklist) {
		this.groundTruthTokenBlacklist = groundTruthTokenBlacklist;
	}

	public String getPathToGoDatabase() {
		return pathToGoDatabase;
	}

	public void setPathToGoDatabase(String pathToGoDatabase) {
		this.pathToGoDatabase = pathToGoDatabase;
	}

	public String getPathToGroundTruthGoAnnotations() {
		return pathToGroundTruthGoAnnotations;
	}

	public void setPathToGroundTruthGoAnnotations(String pathToGroundTruthGoAnnotations) {
		this.pathToGroundTruthGoAnnotations = pathToGroundTruthGoAnnotations;
	}

	public Boolean hasGroundTruthGoAnnotations() {
		return getPathToGroundTruthGoAnnotations() != null && new File(getPathToGroundTruthGoAnnotations()).exists();
	}

	public List<String> getGroundTruthGoAnnotationsFromFile() throws IOException {
		return fromFile(getPathToGroundTruthGoAnnotations());
	}

	public Boolean doCalculateSimpleGoF1Scores() {
		return calculateSimpleGoF1Scores;
	}

	public void setCalculateSimpleGoF1Scores(Boolean calculateSimpleGoF1Scores) {
		this.calculateSimpleGoF1Scores = calculateSimpleGoF1Scores;
	}

	public Boolean doCalculateAncestryGoF1Scores() {
		return calculateAncestryGoF1Scores;
	}

	public void setCalculateAncestryGoF1Scores(Boolean calculateAncestryGoF1Scores) {
		this.calculateAncestryGoF1Scores = calculateAncestryGoF1Scores;
	}

	public Boolean doCalculateSemSimGoF1Scores() {
		return calculateSemSimGoF1Scores;
	}

	public void setCalculateSemSimGoF1Scores(Boolean calculateSemsimGoF1Scores) {
		this.calculateSemSimGoF1Scores = calculateSemsimGoF1Scores;
	}

	public String getPathToGoSlimFile() {
		return pathToGoSlimFile;
	}
	
	public List<String> getGoSlimFile() throws IOException {
		return fromFile(getPathToGoSlimFile());
	}

	public void setPathToGoSlimFile(String pathToGoSlimFile) {
		this.pathToGoSlimFile = pathToGoSlimFile;
	}
	
	public Boolean hasGoSlimFile() {
		return getPathToGoSlimFile() != null && new File(getPathToGoSlimFile()).exists(); 
	}

	public static Pattern getGoSlimFileGotermRegex() {
		return GO_SLIM_FILE_GOTERM_REGEX;
	}

	public int getNumberOfGenerations() {
		return numberOfGenerations;
	}

	public void setNumberOfGenerations(int numberOfGenerations) {
		this.numberOfGenerations = numberOfGenerations;
	}

	public int getPopulationSize() {
		return populationSize;
	}

	public void setPopulationSize(int populationSize) {
		this.populationSize = populationSize;
	}

	public Map<String, Map<String, String>> getCompetitorSettings() {
		return competitorSettings;
	}

	public void setCompetitorSettings(Map<String, Map<String, String>> competitorSettings) {
		this.competitorSettings = competitorSettings;
	}
	
	public Boolean hasCompetitors() {
		return !getCompetitorSettings().isEmpty();
	}
	
	public List<String> getCompetitorDescriptions(String competitor) throws IOException {
		return fromFile(getCompetitorSettings().get(competitor).get(COMPETITOR_DESCRIPTIONS_FILE_KEY));
	}
	
	public List<String> getCompetitorGOAnnotations(String competitor) throws IOException {
		return fromFile(getCompetitorSettings().get(competitor).get(COMPETITOR_GOA_FILE_KEY));
	}

	public boolean doFindHighestPossibleGoScore() {
		return findHighestPossibleGoScore;
	}

	public void setFindHighestPossibleGoScore(boolean findHighestPossibleGoScore) {
		this.findHighestPossibleGoScore = findHighestPossibleGoScore;
	}

	public boolean doWriteFscoreDetailsToOutput() {
		return writeFscoreDetailsToOutput;
	}

	public void setWriteFscoreDetailsToOutput(boolean writeFscoreDetailsToOutput) {
		this.writeFscoreDetailsToOutput = writeFscoreDetailsToOutput;
	}

	public void setFindHighestPossibleEvaluationScore(boolean findHighestPossibleEvaluationScore) {
		this.findHighestPossibleEvaluationScore = findHighestPossibleEvaluationScore;
	}

	public Set<String> getDefaultBlastResultsBlacklist() {
		return defaultBlastResultsBlacklist;
	}

	public void setDefaultBlastResultsBlacklist(Set<String> defaultBlastResultsBlacklist) {
		this.defaultBlastResultsBlacklist = defaultBlastResultsBlacklist;
	}

	public List<String> getDefaultBlastResultsFilter() {
		return defaultBlastResultsFilter;
	}

	public void setDefaultBlastResultsFilter(List<String> defaultBlastResultsFilter) {
		this.defaultBlastResultsFilter = defaultBlastResultsFilter;
	}

	public Set<String> getDefaultTokenBlacklist() {
		return defaultTokenBlacklist;
	}

	public void setDefaultTokenBlacklist(Set<String> defaultTokenBlacklists) {
		this.defaultTokenBlacklist = defaultTokenBlacklists;
	}

	public Pattern getProteinsFastaRegex() {
		return proteinsFastaRegex;
	}

	public void setProteinsFastaRegex(Pattern proteinsFastaRegex) {
		this.proteinsFastaRegex = proteinsFastaRegex;
	}
	
	public Pattern getGroundTruthFastaRegex() {
		return groundTruthFastaRegex;
	}

	public void setGroundTruthFastaRegex(Pattern groundTruthFastaRegex) {
		this.groundTruthFastaRegex = groundTruthFastaRegex;
	}

	public Pattern getSeqSimSearchTableQueryColRegex() {
		return seqSimSearchTableQueryColRegex;
	}

	public void setSeqSimSearchTableQueryColRegex(Pattern seqSimSearchTableQueryColRegex) {
		this.seqSimSearchTableQueryColRegex = seqSimSearchTableQueryColRegex;
	}

	public boolean doFindHighestPossiblePrecision() {
		return findHighestPossiblePrecision;
	}

	public void setFindHighestPossiblePrecision(boolean findHighestPossiblePrecision) {
		this.findHighestPossiblePrecision = findHighestPossiblePrecision;
	}

	public boolean doFindHighestPossibleRecall() {
		return findHighestPossibleRecall;
	}

	public void setFindHighestPossibleRecall(boolean findHighestPossibleRecall) {
		this.findHighestPossibleRecall = findHighestPossibleRecall;
	}

	public boolean doWriteEvaluationSummary() {
		return writeEvaluationSummary;
	}

	public void setWriteEvaluationSummary(boolean writeEvaluationSummary) {
		this.writeEvaluationSummary = writeEvaluationSummary;
	}

	public double getDescriptionScoreThreshold() {
		return this.getParameters().getDescriptionScoreThreshold();
	}

	public void setDescriptionScoreThreshold(double descriptionScoreThreshold) {
		this.getParameters().setDescriptionScoreThreshold(descriptionScoreThreshold);
	}
}