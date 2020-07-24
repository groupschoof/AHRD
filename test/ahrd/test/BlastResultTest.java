package ahrd.test;

import static ahrd.controller.DatabaseSetup.createOrUpdateAhrdDatabase;
import static ahrd.controller.Settings.getSettings;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import org.junit.Before;
import org.junit.Test;

import ahrd.exception.MissingAccessionException;
import ahrd.exception.MissingProteinException;
import ahrd.model.BlastResult;
import ahrd.model.Protein;

public class BlastResultTest {

	public BlastResultTest() {
		super();
	}

	@Before
	public void setUp() throws IOException {
		TestUtils.initTestSettings();
	}

	@Test
	public void testPatternize() {
		BlastResult br = TestUtils.mockBlastResult("accession_5", 5.0, "description_5 Fly-Wing formation", 10, 20, 10,
				20, 200, 30.0, "trembl",
				new HashSet<String>(Arrays.asList("description", "5", "fly", "wing", "formation")));
		assertEquals("5descriptionflyformationwing", br.patternize());
	}

	@Test
	public void testTokenize() throws IOException {
		BlastResult br = new BlastResult("accession_1", 1.0, "one tWo Three protein homolog putative", 10, 20, 10, 20,
				200, 30, "swissprot");
		br.tokenize();
		// test:
		assertEquals(3.0, br.getTokens().size(), 0.0);
		assertTrue(br.getTokens().contains("one"));
		assertTrue(br.getTokens().contains("two"));
		assertTrue(br.getTokens().contains("three"));

		br = new BlastResult("accession_2", 1.0, "Flavohemoprotein-1", 10, 20, 10, 20, 200, 30, "swissprot");
		br.tokenize();
		// test:
		assertEquals(2.0, br.getTokens().size(), 0.0);
		assertTrue(br.getTokens().contains("1"));
		assertTrue(br.getTokens().contains("flavohemoprotein"));
	}

	@Test
	public void testAddBlastResult() throws IOException {
		Map<String, Protein> proteinDb = TestUtils.mockProteinDb();
		Protein p1 = proteinDb.get("gene:chr01.502:mRNA:chr01.502");
		Protein p2 = proteinDb.get("gene:chr01.1056:mRNA:chr01.1056");
		Map<String, List<BlastResult>> blastResults = new HashMap<String, List<BlastResult>>();
		BlastResult.addBlastResult(blastResults,
				new BlastResult("accession_1", 1.0, 10, 20, 10, 20, 200, "swissprot", p1));
		assertEquals(1, blastResults.size());
		assertEquals(1, blastResults.get("accession_1").size());
		BlastResult.addBlastResult(blastResults,
				new BlastResult("accession_1", 1.0, 10, 20, 10, 20, 300, "swissprot", p1));
		assertEquals(1, blastResults.size());
		assertEquals(1, blastResults.get("accession_1").size());
		assertEquals(new Double(300), blastResults.get("accession_1").get(0).getBitScore());
		BlastResult.addBlastResult(blastResults,
				new BlastResult("accession_1", 1.0, 10, 20, 10, 20, 300, "swissprot", p2));
		assertEquals(1, blastResults.size());
		assertEquals(2, blastResults.get("accession_1").size());
		assertEquals(new Double(300), blastResults.get("accession_1").get(1).getBitScore());
		assertEquals(p2, blastResults.get("accession_1").get(1).getProtein());
		BlastResult.addBlastResult(blastResults,
				new BlastResult("accession_2", 1.0, 10, 20, 10, 20, 300, "swissprot", p2));
		assertEquals(2, blastResults.size());
		assertEquals(1, blastResults.get("accession_2").size());
	}

	@Test
	public void testParseBlastResults() throws MissingProteinException, IOException {
		Map<String, Protein> protDb = TestUtils.mockProteinDb();
		Map<String, List<BlastResult>> brs = BlastResult.parseBlastResults(protDb, "tair");
		assertNotNull(brs);
		assertTrue(brs.containsKey("AT3G03300.2"));
		assertTrue(brs.containsKey("AT3G03300.1"));
		assertTrue(brs.containsKey("AT3G43920.1"));
		assertTrue(brs.containsKey("AT3G43920.2"));
		assertTrue(brs.containsKey("AT3G43920.3"));
		assertTrue(brs.containsKey("AT5G20320.1"));
		assertTrue(brs.containsKey("AT1G01040.1"));
		assertEquals(brs.get("AT3G03300.2").get(0).getProtein(), protDb.get("gene:chr01.502:mRNA:chr01.502"));
		assertEquals(brs.get("AT3G03300.1").get(0).getProtein(), protDb.get("gene:chr01.502:mRNA:chr01.502"));
		assertEquals(brs.get("AT3G43920.1").get(0).getProtein(), protDb.get("gene:chr01.502:mRNA:chr01.502"));
		assertEquals(brs.get("AT3G43920.2").get(0).getProtein(), protDb.get("gene:chr01.502:mRNA:chr01.502"));
		assertEquals(brs.get("AT3G43920.3").get(0).getProtein(), protDb.get("gene:chr01.502:mRNA:chr01.502"));
		assertEquals(brs.get("AT5G20320.1").get(0).getProtein(), protDb.get("gene:chr01.502:mRNA:chr01.502"));
		assertEquals(brs.get("AT1G01040.1").get(0).getProtein(), protDb.get("gene:chr01.502:mRNA:chr01.502"));
		BlastResult br = brs.get("AT3G03300.2").get(0);
		assertEquals(br.getAccession(), "AT3G03300.2");
		assertEquals(br.getBitScore(), new Double(94.4), 0.0000001);
		assertEquals(br.getBlastDatabaseName(), "tair");
		assertNull(br.getDescription());
		assertEquals(br.getEValue(), Math.pow(10, -20), Math.pow(10, -21));
		assertEquals(br.getQueryEnd(), new Integer(99));
		assertEquals(br.getQueryStart(), new Integer(1));
		assertEquals(br.getSubjectEnd(), new Integer(1067));
		assertEquals(br.getSubjectStart(), new Integer(969));
		assertNull(br.getSubjectLength());
		// Assert that multiple High Scoring Pairs are read out as a single Hit,
		// that is the one with the best Bit-Score:
		assertEquals(brs.size(), 207);
	}

	@Test
	public void testCompareBlastResultsBasedOnTheirEvalues() {
		BlastResult br1 = new BlastResult("accession_1", 3e-163, "description_1", 10, 20, 10, 20, 200, 30, "swissprot");
		BlastResult br2 = new BlastResult("accession_2", 2e-137, "description_2", 10, 20, 10, 20, 200, 30, "swissprot");
		BlastResult br3 = new BlastResult("accession_3", 2e-137, "description_3", 10, 20, 10, 20, 200, 30, "swissprot");
		assertEquals(-1, br1.compareTo(br2));
		assertEquals(1, br2.compareTo(br1));
		assertEquals(0, br2.compareTo(br3));
	}

	@Test
	public void testFilterBestScoringBlastResults() {
		List<BlastResult> blastResults = new ArrayList<BlastResult>();
		blastResults.add(new BlastResult("accession_1", 3e-163, "description_1", 10, 20, 10, 20, 200, 30, "swissprot"));
		blastResults.add(new BlastResult("accession_2", 3e-163, "description_2", 10, 20, 10, 20, 200, 30, "swissprot"));
		blastResults.add(new BlastResult("accession_3", 2e-137, "description_3", 10, 20, 10, 20, 200, 30, "swissprot"));
		blastResults.add(new BlastResult("accession_4", 2e-137, "description_4", 10, 20, 10, 20, 200, 30, "swissprot"));
		BlastResult br5 = new BlastResult("accession_5", 5.0, "description_5", 10, 20, 10, 20, 200, 30, "swissprot");
		blastResults.add(br5);
		List<BlastResult> fltrdBrs = BlastResult.filterBestScoringBlastResults(blastResults, 4);
		assertNotNull(fltrdBrs);
		assertEquals(4, fltrdBrs.size());
		assertTrue(fltrdBrs.containsAll(blastResults.subList(0, 4)));
		assertTrue(!fltrdBrs.contains(br5));
	}

	@Test
	public void testGetShortAccession() throws IOException {
		BlastResult br = new BlastResult("sp|Q9SXB8|Y1133_ARATH", 1.0, "description_1", 10, 20, 10, 20, 200, 30,
				"swissprot");
		assertEquals("Q9SXB8", br.getShortAccession());
	}

	@Test
	public void testParseLongBlastResults() throws IOException, MissingProteinException {
		TestUtils.initTestSettings();
		getSettings().getBlastDbSettings().get("trembl").put("file", "./test/resources/bgh04634_vs_trEMBL.txt");
		Map<String, Protein> protDb = TestUtils.mockProteinDb();
		Map<String, List<BlastResult>> brs = BlastResult.parseBlastResults(protDb, "trembl");
		assertNotNull(brs);
		assertTrue(!brs.isEmpty());
		// One HSP is double:
		assertEquals(498, brs.size());
		assertEquals(5e-107, brs.get("tr|A0A0G4KNA9|A0A0G4KNA9_9PEZI").get(0).getEValue(), 1e-106);
		assertTrue(brs.containsKey("tr|W9CFB7|W9CFB7_9HELO"));
		// Verify the BlastResult is valid:
		assertEquals(1, brs.get("tr|W9CFB7|W9CFB7_9HELO").size());
		BlastResult br = brs.get("tr|W9CFB7|W9CFB7_9HELO").get(0);
		assertNotNull(br);
		assertNotNull(br.getAccession());
		assertNotNull(br.getBitScore());
		assertNotNull(br.getQueryEnd());
		assertNotNull(br.getQueryStart());
		assertNotNull(br.getEValue());
		assertEquals(Math.pow(10, -163), br.getEValue(), Math.pow(10, -162));
		assertNotNull(br.getBlastDatabaseName());
		assertNotNull(br.getSubjectEnd());
		assertNotNull(br.getSubjectStart());
		// The following fields are not yet set, because they have to be read
		// out of the FASTA database:
		assertNull(br.getDescription());
		assertNull(br.getSubjectLength());
		assertNotNull(br.getTokens());
		assertTrue(br.getTokens().isEmpty());
		// Because of the missing field values the BlastResult should not yet be
		// valid:
		assertTrue(!br.isValid());
	}

	@Test
	public void testParseLongBlastDatabase() throws IOException, MissingProteinException, MissingAccessionException {
		TestUtils.initTestSettings();
		getSettings().getBlastDbSettings().get("trembl").put("file", "./test/resources/bgh04634_vs_trEMBL.txt");
		getSettings().getBlastDbSettings().get("trembl").put("database",
				"./test/resources/bgh04634_trembl_database.fasta");
		getSettings().setPathToGeneOntologyResults(null);
		createOrUpdateAhrdDatabase(false);
		Map<String, Protein> protDb = TestUtils.mockProteinDb();
		Map<String, List<BlastResult>> brs = BlastResult.parseBlastResults(protDb, "trembl");
		BlastResult.parseBlastDatabase(brs);
		Protein p1 = protDb.get("gene:chr01.1056:mRNA:chr01.1056");
		assertTrue(!p1.getBlastResults().get("trembl").isEmpty());
		Map<String, BlastResult> trBrs = new HashMap<String, BlastResult>();
		for (BlastResult trBr : p1.getBlastResults().get("trembl")) {
			trBrs.put(trBr.getAccession(), trBr);
		}
		// System.out.println(trAccs);
		assertTrue(trBrs.containsKey("tr|W9CFB7|W9CFB7_9HELO"));
		BlastResult theBr = trBrs.get("tr|W9CFB7|W9CFB7_9HELO");
		assertNotNull(theBr);
		assertNotNull(theBr.getSubjectLength());
		assertNotNull(theBr.getDescription());
		assertNotNull(theBr.getTokens());
		assertEquals(3, theBr.getTokens().size());
		assertTrue(theBr.isValid());
	}

	@Test
	public void testFilterBestScoringBlastResultsLong()
			throws IOException, MissingProteinException, MissingAccessionException {
		TestUtils.initTestSettings();
		getSettings().getBlastDbSettings().get("trembl").put("file", "./test/resources/bgh04634_vs_trEMBL.txt");
		getSettings().getBlastDbSettings().get("trembl").put("database",
				"./test/resources/bgh04634_trembl_database.fasta");
		getSettings().setPathToGeneOntologyResults(null);
		Map<String, Protein> protDb = TestUtils.mockProteinDb();
		Map<String, List<BlastResult>> brs = BlastResult.parseBlastResults(protDb, "trembl");
		BlastResult.parseBlastDatabase(brs);
		Protein p1 = protDb.get("gene:chr01.1056:mRNA:chr01.1056");
		assertTrue(!p1.getBlastResults().get("trembl").isEmpty());
		List<BlastResult> filtrdBrs = BlastResult.filterBestScoringBlastResults(p1.getBlastResults().get("trembl"),
				200);
		assertNotNull(filtrdBrs);
		assertEquals(200, filtrdBrs.size());
		assertEquals("tr|W9CFB7|W9CFB7_9HELO", filtrdBrs.get(0).getAccession());
	}

	@Test
	public void testGenerateHRDCandidateForProtein()
			throws IOException, MissingProteinException, MissingAccessionException {
		TestUtils.initTestSettings();
		getSettings().getBlastDbSettings().get("trembl").put("file", "./test/resources/bgh04634_vs_trEMBL.txt");
		getSettings().getBlastDbSettings().get("trembl").put("database",
				"./test/resources/bgh04634_trembl_database.fasta");
		getSettings().setPathToGeneOntologyResults(null);
		getSettings().setWriteBestBlastHitsToOutput(true);
		Map<String, Protein> protDb = TestUtils.mockProteinDb();
		Map<String, List<BlastResult>> brs = BlastResult.parseBlastResults(protDb, "trembl");
		BlastResult.parseBlastDatabase(brs);
		Protein p1 = protDb.get("gene:chr01.1056:mRNA:chr01.1056");
		assertTrue(!p1.getEvaluationScoreCalculator().getUnchangedBlastResults().isEmpty());
		assertNotNull(p1.getEvaluationScoreCalculator().getUnchangedBlastResults().get("trembl"));
		assertEquals("tr|W9CFB7|W9CFB7_9HELO",
				p1.getEvaluationScoreCalculator().getUnchangedBlastResults().get("trembl").getAccession());
	}
}