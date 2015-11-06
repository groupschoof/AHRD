package ahrd.test;

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

import org.junit.Test;

import ahrd.exception.MissingProteinException;
import ahrd.model.BlastResult;
import ahrd.model.Protein;

public class BlastResultTest {

	public BlastResultTest() {
		super();
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
		TestUtils.initTestSettings();
		BlastResult br = new BlastResult("accession_1", 1.0, "one tWo Three protein homolog putative", 10, 20, 10, 20,
				200, 30, "swissprot");
		br.tokenize();
		// test:
		assertEquals(3.0, br.getTokens().size(), 0.0);
		assertTrue(br.getTokens().contains("one"));
		assertTrue(br.getTokens().contains("two"));
		assertTrue(br.getTokens().contains("three"));
	}

	@Test
	public void testAddBlastResult() throws IOException {
		TestUtils.initTestSettings();
		Map<String, Protein> proteinDb = TestUtils.mockProteinDb();
		Protein p1 = proteinDb.get("gene:chr01.502:mRNA:chr01.502");
		Protein p2 = proteinDb.get("gene:chr01.1056:mRNA:chr01.1056");
		Map<String, List<BlastResult>> blastResults = new HashMap<String, List<BlastResult>>();
		BlastResult.addBlastResult(blastResults,
				new BlastResult("accession_1", 1.0, 10, 20, 10, 20, 200, "swissprot", p1), null);
		assertEquals(1, blastResults.size());
		assertEquals(1, blastResults.get("accession_1").size());
		BlastResult.addBlastResult(blastResults,
				new BlastResult("accession_1", 1.0, 10, 20, 10, 20, 300, "swissprot", p1), null);
		assertEquals(1, blastResults.size());
		assertEquals(1, blastResults.get("accession_1").size());
		assertEquals(new Double(300), blastResults.get("accession_1").get(0).getBitScore());
		BlastResult.addBlastResult(blastResults,
				new BlastResult("accession_1", 1.0, 10, 20, 10, 20, 300, "swissprot", p2), null);
		assertEquals(1, blastResults.size());
		assertEquals(2, blastResults.get("accession_1").size());
		assertEquals(new Double(300), blastResults.get("accession_1").get(1).getBitScore());
		assertEquals(p2, blastResults.get("accession_1").get(1).getProtein());
		BlastResult.addBlastResult(blastResults,
				new BlastResult("accession_2", 1.0, 10, 20, 10, 20, 300, "swissprot", p2), null);
		assertEquals(2, blastResults.size());
		assertEquals(1, blastResults.get("accession_2").size());
	}

	@Test
	public void testParseBlastResults() throws MissingProteinException, IOException {
		TestUtils.initTestSettings();
		Map<String, Protein> protDb = TestUtils.mockProteinDb();
		Map<String, List<BlastResult>> brs = BlastResult.parseBlastResults(protDb, "tair", null);
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
	public void testParseBlastDatabase() throws IOException, MissingProteinException {
		TestUtils.initTestSettings();
		Map<String, Protein> protDb = TestUtils.mockProteinDb();
		Map<String, List<BlastResult>> brs = BlastResult.parseBlastResults(protDb, "tair", null);
		BlastResult.parseBlastDatabase(protDb, "tair", brs);
		Protein p1 = protDb.get("gene:chr01.502:mRNA:chr01.502");
		Protein p2 = protDb.get("gene:chr01.1056:mRNA:chr01.1056");
		assertTrue(!p1.getBlastResults().get("tair").isEmpty());
		assertEquals(7, p1.getBlastResults().get("tair").size());
		assertEquals("AT3G03300.2", p1.getBlastResults().get("tair").get(0).getAccession());
		assertEquals(new Integer(1375), p1.getBlastResults().get("tair").get(0).getSubjectLength());
		assertTrue(!p2.getBlastResults().get("tair").isEmpty());
		assertEquals(200, p2.getBlastResults().get("tair").size());
		assertEquals("AT3G45420.1", p2.getBlastResults().get("tair").get(199).getAccession());
		assertEquals(new Integer(668), p2.getBlastResults().get("tair").get(199).getSubjectLength());
	}

	@Test
	public void testCompareBlastResultsBasedOnTheirEvalues() {
		BlastResult br1 = new BlastResult("accession_1", 1.0, "description_1", 10, 20, 10, 20, 200, 30, "swissprot");
		BlastResult br2 = new BlastResult("accession_2", 2.0, "description_2", 10, 20, 10, 20, 200, 30, "swissprot");
		BlastResult br3 = new BlastResult("accession_3", 2.0, "description_3", 10, 20, 10, 20, 200, 30, "swissprot");
		assertEquals(-1, br1.compareTo(br2));
		assertEquals(1, br2.compareTo(br1));
		assertEquals(0, br2.compareTo(br3));
	}

	@Test
	public void testFilterBestScoringBlastResults() {
		List<BlastResult> blastResults = new ArrayList<BlastResult>();
		blastResults.add(new BlastResult("accession_1", 1.0, "description_1", 10, 20, 10, 20, 200, 30, "swissprot"));
		blastResults.add(new BlastResult("accession_2", 2.0, "description_2", 10, 20, 10, 20, 200, 30, "swissprot"));
		blastResults.add(new BlastResult("accession_3", 3.0, "description_3", 10, 20, 10, 20, 200, 30, "swissprot"));
		BlastResult br = new BlastResult("accession_4", 4.0, "description_4", 10, 20, 10, 20, 200, 30, "swissprot");
		blastResults.add(br);
		blastResults.add(new BlastResult("accession_5", 5.0, "description_5", 10, 20, 10, 20, 200, 30, "swissprot"));
		blastResults = BlastResult.filterBestScoringBlastResults(blastResults, 4);
		assertEquals(4, blastResults.size());
		assertTrue(blastResults.contains(br));
	}

	@Test
	public void testGetShortAccession() {
		BlastResult br = new BlastResult("sp|Q9SXB8|Y1133_ARATH", 1.0, "description_1", 10, 20, 10, 20, 200, 30,
				"swissprot");
		assertEquals("Q9SXB8", br.getShortAccession());
	}

	@Test
	public void testParseLongBlastResults() throws IOException, MissingProteinException {
		TestUtils.initTestSettings();
		getSettings().getBlastDbSettings().get("trembl").put("file", "./test/resources/bgh04634_vs_trEMBL.txt");
		Map<String, Protein> protDb = TestUtils.mockProteinDb();
		Map<String, List<BlastResult>> brs = BlastResult.parseBlastResults(protDb, "trembl", null);
		assertNotNull(brs);
		assertTrue(brs.containsKey("tr|W9CFB7|W9CFB7_9HELO"));
	}

	@Test
	public void testParseLongBlastDatabase() throws IOException, MissingProteinException {
		TestUtils.initTestSettings();
		getSettings().getBlastDbSettings().get("trembl").put("file", "./test/resources/bgh04634_vs_trEMBL.txt");
		getSettings().getBlastDbSettings().get("trembl").put("database",
				"/projects/irg/grp_pcb/projects/AHRD/blast_dbs_15_Oct_2015/uniprot_trembl.fasta");
		Map<String, Protein> protDb = TestUtils.mockProteinDb();
		Map<String, List<BlastResult>> brs = BlastResult.parseBlastResults(protDb, "trembl", null);
		BlastResult.parseBlastDatabase(protDb, "trembl", brs);
		Protein p1 = protDb.get("gene:chr01.1056:mRNA:chr01.1056");
		assertTrue(!p1.getBlastResults().get("tair").isEmpty());
		assertEquals(200, p1.getBlastResults().get("trembl").size());
		assertEquals("tr|W9CFB7|W9CFB7_9HELO", p1.getBlastResults().get("trembl").get(0).getAccession());
	}
}
