package ahrd.test;

import static ahrd.controller.Parameters.randomParameters;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.junit.Before;
import org.junit.Test;

import ahrd.controller.Batcher;
import ahrd.controller.Settings;
import ahrd.controller.TrainerBatcher;

import com.esotericsoftware.yamlbeans.YamlException;
import com.esotericsoftware.yamlbeans.YamlReader;

public class TrainerBatcherTest {

	private TrainerBatcher trainerBatcher;

	public TrainerBatcherTest() {
		super();
	}

	@Before
	public void setUp() throws IOException {
		this.trainerBatcher = new TrainerBatcher(mockInput());
	}

	@SuppressWarnings("unchecked")
	private Map<String, Object> mockInput() throws FileNotFoundException,
			YamlException {
		YamlReader reader = new YamlReader(new FileReader(
				"./test/resources/trainer_batcher_input_test.yml"));
		return (Map<String, Object>) reader.read();
	}

	@Test
	public void testBatch() {
		this.trainerBatcher.batch();
		assertEquals(10, this.trainerBatcher.getOutput().size());
	}

	@Test
	@SuppressWarnings("unchecked")
	public void testGenerateYml() {
		List<String> blastDbs = new ArrayList<String>();
		blastDbs.add("swissprot");
		blastDbs.add("tair");
		blastDbs.add("trembl");
		Map<String, Object> batchYml = this.trainerBatcher.generateYml(
				"testSimAnnealBatch.yml", randomParameters(blastDbs));
		// batchYml well formed?
		assertTrue(batchYml != null);
		assertTrue("Generated Yml-Hash should contain blast_dbs.",
				batchYml.containsKey(Settings.BLAST_DBS_KEY));
		// Isn't type-casting a nice way to maintain yourself busy?!
		Map<String, String> sprotBlastDb = (Map<String, String>) ((Map<String, Object>) batchYml
				.get(Settings.BLAST_DBS_KEY)).get("swissprot");
		assertTrue(
				"Blast-Database-Weight of swissprot is null or empty String!",
				assertNotNullAndNotEmpty(sprotBlastDb
						.get(Settings.BLAST_DB_WEIGHT_KEY)));
		assertEquals("./test/resources/blacklist_descline.txt",
				sprotBlastDb.get(Settings.BLAST_BLACKLIST_KEY));
		assertEquals("./test/resources/filter_descline_sprot.txt",
				sprotBlastDb.get(Settings.BLAST_FILTER_KEY));
		assertEquals("./test/resources/blacklist_token.txt",
				sprotBlastDb.get(Settings.TOKEN_BLACKLIST_KEY));
		assertEquals("./test/resources/swissprot.pairwise",
				sprotBlastDb.get(Settings.BLAST_RESULT_FILE_KEY));
		// Interpro:
		assertEquals("./test/resources/interpro_31.xml",
				batchYml.get(Settings.INTERPRO_DATABASE_KEY).toString());
		assertEquals("./test/resources/interpro_result.raw",
				batchYml.get(Settings.INTERPRO_RESULT_KEY).toString());
		// Gene-Ontology:
		assertEquals("./test/resources/go_results.csv",
				batchYml.get(Settings.GENE_ONTOLOGY_RESULT_KEY).toString());
		// Test Output-File:
		assertEquals(
				"./test/resources/testSimAnnealBatch_ahrd_trainer_out.csv",
				batchYml.get(Settings.OUTPUT_KEY));
		// Test existence of 'Path to Batch-Yml':
		// - Note: The correctness of the path is tested elsewhere.
		assertTrue("Batch's YML should hold the path to itself.",
				batchYml.containsKey(Batcher.PATH_TO_BATCH_YML_KEY));
		// Assert proteins-fasta:
		assertEquals("./test/resources/proteins.fasta",
				batchYml.get(Settings.PROTEINS_FASTA_KEY));
		// Assert remember_simulated_annealing_path:
		assertTrue(
				"remember_simulated_annealing_path should not be null",
				batchYml.get(Settings.REMEMBER_SIMULATED_ANNEALING_PATH_KEY) != null);
		assertTrue("remember_simulated_annealing_path should be true",
				Boolean.parseBoolean(batchYml.get(
						Settings.REMEMBER_SIMULATED_ANNEALING_PATH_KEY)
						.toString()));
		// Assert boolean parameter 'find_highest_possible_evaluation_score'
		// gets passed on:
		assertNotNull(
				"Batch's YML should contain boolean parameter 'find_highest_possible_evaluation_score'",
				batchYml.get(Settings.FIND_HIGHEST_POSSIBLE_EVALUATION_SCORE_KEY));
		assertEquals(
				"Batch's YML should have set boolean parameter 'find_highest_possible_evaluation_score' to TRUE.",
				"true",
				batchYml.get(
						Settings.FIND_HIGHEST_POSSIBLE_EVALUATION_SCORE_KEY)
						.toString());
	}

	private boolean assertNotNullAndNotEmpty(String toBeValidated) {
		return toBeValidated != null && !toBeValidated.equals("");
	}
}
