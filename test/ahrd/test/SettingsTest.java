package ahrd.test;

import static ahrd.controller.Settings.getSettings;
import static ahrd.controller.Settings.setSettings;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.util.regex.Pattern;

import org.junit.Before;
import org.junit.Test;

import ahrd.controller.Settings;

public class SettingsTest {

	@Before
	public void setUp() throws IOException {
		TestUtils.initTestSettings();
	}

	@Test
	public void testAhrdLoadsYamlInput() throws IOException {
		assertTrue(getSettings() instanceof Settings);
		assertNotNull("Input-Map should have proteins_fasta", getSettings()
				.getProteinsFasta());
		assertEquals(3.0, getSettings().getBlastDatabases().size(), 0.0);
		assertNotNull(getSettings().getBlastDatabases().contains("swissprot"));
		assertNotNull(getSettings().getBlastDatabases().contains("tair"));
		assertNotNull(getSettings().getBlastDatabases().contains("trembl"));
		assertNotNull(getSettings().getTokenBlacklist("swissprot"));
		assertNotNull(getSettings().getTokenBlacklist("tair"));
		assertNotNull(getSettings().getTokenBlacklist("trembl"));
		assertEquals(
				"Test-Temperature in input.yml is set to 10 and should have been initialized correctly.",
				Integer.valueOf(10), getSettings().getTemperature());
	}

	@Test
	public void testAhrdLoadsSequenceSimilaritySearchTableSettings()
			throws IOException {
		// Assert default values:
		assertNull(getSettings().getSeqSimSearchTableCommentLineRegex());
		assertEquals("\t", getSettings().getSeqSimSearchTableSep());
		assertEquals(Integer.valueOf(0), getSettings()
				.getSeqSimSearchTableQueryCol());
		assertEquals(Integer.valueOf(1), getSettings()
				.getSeqSimSearchTableSubjectCol());
		assertEquals(Integer.valueOf(6), getSettings()
				.getSeqSimSearchTableQueryStartCol());
		assertEquals(Integer.valueOf(7), getSettings()
				.getSeqSimSearchTableQueryEndCol());
		assertEquals(Integer.valueOf(8), getSettings()
				.getSeqSimSearchTableSubjectStartCol());
		assertEquals(Integer.valueOf(9), getSettings()
				.getSeqSimSearchTableSubjectEndCol());
		assertEquals(Integer.valueOf(10), getSettings()
				.getSeqSimSearchTableEValueCol());
		assertEquals(Integer.valueOf(11), getSettings()
				.getSeqSimSearchTableBitScoreCol());
		assertEquals(
				Pattern.compile(
						"^>(?<accession>[aA][tT][0-9mMcC][gG]\\d+(\\.\\d+)?)\\s+\\|[^\\|]+\\|\\s+(?<description>[^\\|]+)(\\s*\\|.*)?$")
						.toString(), getSettings().getFastaHeaderRegex("tair")
						.toString());
		assertEquals(Settings.DEFAULT_FASTA_HEADER_REGEX.toString(),
				getSettings().getFastaHeaderRegex("trembl").toString());
		// Assert custom values:
		setSettings(new Settings(
				"./test/resources/ahrd_input_seq_sim_table.yml"));
		assertEquals(Pattern.compile("#").toString(), getSettings()
				.getSeqSimSearchTableCommentLineRegex().toString());
		assertEquals("\t", getSettings().getSeqSimSearchTableSep());
		assertEquals(Integer.valueOf(10), getSettings()
				.getSeqSimSearchTableQueryCol());
		assertEquals(Integer.valueOf(11), getSettings()
				.getSeqSimSearchTableSubjectCol());
		assertEquals(Integer.valueOf(16), getSettings()
				.getSeqSimSearchTableQueryStartCol());
		assertEquals(Integer.valueOf(17), getSettings()
				.getSeqSimSearchTableQueryEndCol());
		assertEquals(Integer.valueOf(18), getSettings()
				.getSeqSimSearchTableSubjectStartCol());
		assertEquals(Integer.valueOf(19), getSettings()
				.getSeqSimSearchTableSubjectEndCol());
		assertEquals(Integer.valueOf(20), getSettings()
				.getSeqSimSearchTableEValueCol());
		assertEquals(Integer.valueOf(21), getSettings()
				.getSeqSimSearchTableBitScoreCol());
		assertEquals(
				Pattern.compile(
						"^>(?<accession>[aA][tT][0-9mMcC][gG]\\d+(\\.\\d+)?)\\s+\\|[^\\|]+\\|\\s+(?<description>[^\\|]+)(\\s*\\|.*)?$")
						.toString(), getSettings().getFastaHeaderRegex("tair")
						.toString());
		assertEquals(Settings.DEFAULT_FASTA_HEADER_REGEX.toString(),
				getSettings().getFastaHeaderRegex("trembl").toString());
	}

	@Test
	public void testClone() {
		Settings s = getSettings();
		s.setAvgDescriptionEvaluationScore(0.8);
		Settings c = s.clone();
		// test
		assertTrue("A clone should not be it's 'parent'.", s != c);
		// All primitive types get automatically cloned by Object.clone(),
		// so just test the Maps:
		for (String blastDb : s.getBlastDatabases()) {
			assertTrue("The BlastDatabaseMaps should not be the same objects.",
					s.getBlastDbSettings().get(blastDb) != c.getBlastDbSettings().get(blastDb));
			assertTrue("The BlastDatabaseMaps should not have the same identityHashCode.",
					System.identityHashCode(s.getBlastDbSettings().get(blastDb))
					!= System.identityHashCode(c.getBlastDbSettings().get(blastDb)));
		}
		// Test passing on the average evaluation score:
		assertEquals(s.getAvgDescriptionEvaluationScore(), c.getAvgDescriptionEvaluationScore(), 0.0);
		// Assure they are different Objects. As the operator != does not reveal
		// this, we set the clone to a different value than its parent:
		c.setAvgDescriptionEvaluationScore(0.7);
		assertTrue(
				"Cloning should result in the average evaluation scores being different Objects.",
				s.getAvgDescriptionEvaluationScore() != c.getAvgDescriptionEvaluationScore());
	}

	@Test
	public void testRememberSimulatedAnnealingPath() {
		// Flag should be set in default Test-Settings:
		assertTrue(
				"Path through Parameter-Space should be remembered, but flag is set to FALSE.",
				getSettings().rememberSimulatedAnnealingPath());
	}
	
	@Test
	public void testEvaluateOnlyValidTokens() {
		assertTrue(
				"The flag evaluateOnlyValidTokens should be set to TRUE by default, but it is set to FALSE.",
				getSettings().getEvaluateOnlyValidTokens());
	}
	
}
