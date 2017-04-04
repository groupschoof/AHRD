package ahrd.test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.assertNotNull;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.Map;

import org.junit.Before;
import org.junit.Test;

import ahrd.controller.Batcher;
import ahrd.controller.Settings;

import com.esotericsoftware.yamlbeans.YamlException;
import com.esotericsoftware.yamlbeans.YamlReader;

public class BatcherTest {

	private Batcher batcher;

	public BatcherTest() {
		super();
	}

	@Before
	public void setUp() throws FileNotFoundException, YamlException {
		this.batcher = new Batcher(mockInput());
	}

	@SuppressWarnings("unchecked")
	private Map<String, Object> mockInput() throws FileNotFoundException,
			YamlException {
		YamlReader reader = new YamlReader(new FileReader(
				"./test/resources/batcher_input_test.yml"));
		return (Map<String, Object>) reader.read();
	}

	@Test
	public void testFindFileInDirectory() {
		String batchName = "go_references.fasta";
		String fileToBeFound = "./test/resources/go_references.csv";
		String directory = "./test/resources/";
		assertEquals(fileToBeFound,
				this.batcher.findFileInDirectory(directory, batchName));
	}

	@Test
	public void testBatch() {
		this.batcher.batch();
		assertTrue(
				"After batching Batcher should hold HashMaps representing the AHRD-Input-YMLs.",
				!this.batcher.getOutput().isEmpty());
		assertEquals(3, this.batcher.getOutput().size());
	}

	@Test
	@SuppressWarnings("unchecked")
	public void testGenerateYml() throws FileNotFoundException, YamlException {
		Map<String, Object> batchYml = this.batcher
				.generateYml("batch001.fasta");

		assertTrue("Generated Yml-Hash should contain blast_dbs.",
				batchYml.containsKey(Settings.BLAST_DBS_KEY));
		// Isn't type-casting a nice way to maintain yourself busy?!
		Map<String, String> sprotBlastDb = (Map<String, String>) ((Map<String, Object>) batchYml
				.get(Settings.BLAST_DBS_KEY)).get("swissprot");
		assertEquals("100", sprotBlastDb.get(Settings.BLAST_DB_WEIGHT_KEY));
		assertEquals("./test/resources/filter_descline_sprot.txt",
				sprotBlastDb.get(Settings.BLAST_FILTER_KEY));
		assertEquals("./test/resources/sprot_blast_results/batch001.pairwise",
				sprotBlastDb.get(Settings.BLAST_RESULT_FILE_KEY));
		assertEquals("./test/resources/swissprot_blast_db.fasta",
				sprotBlastDb.get(Settings.BLAST_DATABASE_KEY));
		// Gene-Ontology:
		assertEquals("./test/resources/gene_ontology_references/batch001.csv",
				sprotBlastDb.get(Settings.GENE_ONTOLOGY_REFERENCE_KEY).toString());
		// Verify, that optional FASTA_HEADER_REGEX parameters are also passed
		// on:
		assertEquals(
				"^>(?<accession>[aA][tT][0-9mMcC][gG]\\d+(\\.\\d+)?)\\s+\\|[^\\|]+\\|\\s+(?<description>[^\\|]+)(\\s*\\|.*)?$",
				((Map<String, String>) ((Map<String, Object>) batchYml
						.get(Settings.BLAST_DBS_KEY)).get("tair"))
						.get(Settings.FASTA_HEADER_REGEX_KEY));
		// Interpro:
		assertEquals("./test/resources/interpro_31.xml",
				batchYml.get(Settings.INTERPRO_DATABASE_KEY).toString());
		assertEquals("./test/resources/interpro_results/batch001.raw", batchYml
				.get(Settings.INTERPRO_RESULT_KEY).toString());
		// Test Output-File:
		assertEquals("./test/resources/batch001_ahrd_out.csv",
				batchYml.get(Settings.OUTPUT_KEY));
		// Test existence of 'Path to Batch-Yml':
		// - Note: The correctness of the path is tested elsewhere.
		assertTrue("Batch's YML should hold the path to itself.",
				batchYml.containsKey(Batcher.PATH_TO_BATCH_YML_KEY));
		// Assert proteins-fasta:
		assertEquals("./test/resources/proteins/batch001.fasta",
				batchYml.get(Settings.PROTEINS_FASTA_KEY));
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

	@Test
	public void testGeneratePathToBatchYml() {
		String batchName = "oveja001.fasta";
		assertEquals("./test/resources/batch_ymls/oveja001.yml",
				this.batcher.generatePathToBatchYml(batchName));
	}
}