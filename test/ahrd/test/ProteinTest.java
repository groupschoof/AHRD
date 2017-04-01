package ahrd.test;

import static ahrd.controller.Settings.getSettings;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

import org.junit.Before;
import org.junit.Test;

import ahrd.exception.MissingAccessionException;
import ahrd.model.Protein;

public class ProteinTest {

	public ProteinTest() {
		super();
	}

	@Before
	public void setup() throws IOException {
		TestUtils.initTestSettings();
	}

	@Test
	public void testProtConstrInitsTrainScorCalcWhenNeeded() throws IOException {
		TestUtils.initTestSettings();
		Protein p = TestUtils.mockProtein();
		// The above Test-Settings are set to Evaluator-Mode and set to Write
		// out
		// the best Blast-Hits, so we should have a non-null
		// EvaluationScoreCalculator:
		assertTrue(
				"Expected mocked protein's EvaluationScoreCalculator not te be null!",
				p.getEvaluationScoreCalculator() != null);
		// Test the opposite:
		getSettings().setWriteBestBlastHitsToOutput(false);
		getSettings().setPathToGroundTruthFasta(null);
		p = TestUtils.mockProtein();
		assertTrue(
				"Expected mocked protein's EvaluationScoreCalculator te be null!",
				p.getEvaluationScoreCalculator() == null);
	}

	@Test
	public void testProteinConstructionFromFasta()
			throws MissingAccessionException {
		// NOTE: We expect this method to be called with a split String, read
		// from a FASTA-File.
		// The splitting is done for '>'!
		String input = "MySequence | more information | we don't need\n"
				+ "MADDSKFCFFLVSTFLLLAVVVNVTLAANYVPGDDILLNCGGPDNLPDADGRKWGTDIGS\n"
				+ "KYMLGSKSSTSDAADQKPSVPQVPFMSARVFQSEFTYSFPVASGRKFVRLYFYPSSYNKL\n"
				+ "NATNAIFNVASGQYTLLRNFSAAQTAEALNYDYLTKEFSINVRWKYLHCNVISPSGDTGM\n"
				+ "SPGYDASMTDSRSSGISMSIGGRSLASEDSDGLTPSAVFSQIMNPKGR*";
		String sequence = "MADDSKFCFFLVSTFLLLAVVVNVTLAANYVPGDDILLNCGGPDNLPDADGRKWGTDIGS"
				+ "KYMLGSKSSTSDAADQKPSVPQVPFMSARVFQSEFTYSFPVASGRKFVRLYFYPSSYNKL"
				+ "NATNAIFNVASGQYTLLRNFSAAQTAEALNYDYLTKEFSINVRWKYLHCNVISPSGDTGM"
				+ "SPGYDASMTDSRSSGISMSIGGRSLASEDSDGLTPSAVFSQIMNPKGR*";
		Protein prot = Protein.constructFromFastaEntry(input);
		assertNotNull(
				"Protein should be constructable from a single FASTA-Entry.",
				prot);
		assertEquals(
				"Constructed protein should have it's Accession correctly set",
				"MySequence", prot.getAccession());
		assertEquals(((new Integer(sequence.length())).longValue()),
				(prot.getSequenceLength().longValue()));
		// In order to save memory, when not writing output in fasta-format,
		// ensure the AA-sequences are not held in memory:
		assertNull(
				"The AA-Sequence should not be held in memory, as output-format is NOT set to fasta.",
				prot.getSequence());
	}

	@Test
	public void testProteinConstructionFromFastaWithSequence()
			throws MissingAccessionException, IOException {
		// Setup trigger to also remember the AA-Sequence in memory:
		getSettings().setOutputFasta(true);
		// NOTE: We expect this method to be called with a split String, read
		// from a FASTA-File.
		// The splitting is done for '>'!
		String input = "MySequence | more information | we don't need\n"
				+ "MADDSKFCFFLVSTFLLLAVVVNVTLAANYVPGDDILLNCGGPDNLPDADGRKWGTDIGS\n"
				+ "KYMLGSKSSTSDAADQKPSVPQVPFMSARVFQSEFTYSFPVASGRKFVRLYFYPSSYNKL\n"
				+ "NATNAIFNVASGQYTLLRNFSAAQTAEALNYDYLTKEFSINVRWKYLHCNVISPSGDTGM\n"
				+ "SPGYDASMTDSRSSGISMSIGGRSLASEDSDGLTPSAVFSQIMNPKGR*";
		String sequence = "MADDSKFCFFLVSTFLLLAVVVNVTLAANYVPGDDILLNCGGPDNLPDADGRKWGTDIGS"
				+ "KYMLGSKSSTSDAADQKPSVPQVPFMSARVFQSEFTYSFPVASGRKFVRLYFYPSSYNKL"
				+ "NATNAIFNVASGQYTLLRNFSAAQTAEALNYDYLTKEFSINVRWKYLHCNVISPSGDTGM"
				+ "SPGYDASMTDSRSSGISMSIGGRSLASEDSDGLTPSAVFSQIMNPKGR*";
		Protein prot = Protein.constructFromFastaEntry(input);
		assertNotNull(
				"Protein should be constructable from a single FASTA-Entry.",
				prot);
		assertEquals(
				"Constructed protein should have it's Accession correctly set",
				"MySequence", prot.getAccession());
		assertEquals(((new Integer(sequence.length())).longValue()),
				(prot.getSequenceLength().longValue()));
		assertNotNull("Protein should have its AA-sequence set.",
				prot.getSequence());
		assertEquals(
				"The Protein's AA-Sequence mustn't be changed, except for stripped newlines.",
				sequence, prot.getSequence());
	}

	@Test
	public void testInitialisationOfProteinsFromFasta() throws IOException,
			MissingAccessionException {
		TestUtils.initTestSettings();
		Map<String, Protein> prot_db = null;
		prot_db = Protein.initializeProteins(getSettings().getProteinsFasta());
		assertTrue(prot_db instanceof Map);
		assertNotNull(prot_db);
		assertEquals(3, prot_db.size());
	}

	@Test
	public void testFastaSplitter() {
		String str = ">Protein One->Expressed in stupid user\nSPGYDASMTDSRSSGISMSIGGRSLASEDSDGLTPSAVFSQIMNPKGR\n>Protein Two\nMADDSKFCFFLVSTFLLLAVVVNVTLAANYVPGDDILLNCGGPDNLPDADGRKWGTDIGS";
		List<String> fastaEntries = Protein.splitFasta(str);
		assertEquals(
				Arrays.asList(
						"Protein One->Expressed in stupid user\nSPGYDASMTDSRSSGISMSIGGRSLASEDSDGLTPSAVFSQIMNPKGR",
						"Protein Two\nMADDSKFCFFLVSTFLLLAVVVNVTLAANYVPGDDILLNCGGPDNLPDADGRKWGTDIGS"),
				fastaEntries);
	}

}
