package ahrd.test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.junit.Before;
import org.junit.Test;
import org.xml.sax.SAXException;

import ahrd.controller.AHRD;
import ahrd.exception.MissingAccessionException;
import ahrd.exception.MissingProteinException;
import ahrd.model.BlastResult;
import ahrd.model.DomainScoreCalculator;
import ahrd.model.InterproResult;
import ahrd.model.Protein;

public class AhrdTest {

	private AHRD ahrd;

	public AhrdTest() {
	}

	@Before
	public void setUp() throws IOException {
		ahrd = new AHRD("./test/resources/ahrd_input.yml");
	}

	@Test
	public void testAhrdInitializesProteins() throws IOException,
			MissingAccessionException {
		ahrd.initializeProteins();
		assertNotNull(ahrd.getProteins());
		assertTrue(ahrd.getProteins().containsKey(
				"gene:chr01.502:mRNA:chr01.502"));
		assertTrue(ahrd.getProteins().containsKey(
				"gene:chr01.1056:mRNA:chr01.1056"));
	}

	@Test
	public void testAhrdParsesBlast() throws IOException,
			MissingProteinException, SAXException, MissingAccessionException {
		// Start Setup the test!
		// We need the test-protein-Database in memory:
		ahrd.setProteins(TestUtils.mockProteinDb());
		// now, we can parse the BlastResults:
		ahrd.parseBlastResults();
		for (String iterProtAcc : ahrd.getProteins().keySet()) {
			Protein protein = ahrd.getProteins().get(iterProtAcc);
			assertTrue(protein.getBlastResults().size() > 0);
			assertTrue(protein.getBlastResults().containsKey("swissprot"));
			assertTrue(protein.getBlastResults().containsKey("tair"));
			assertTrue(protein.getBlastResults().containsKey("trembl"));
		}
		// Start the actual test:
		// test number of tokens:
		Protein p = ahrd.getProteins().get("gene:chr01.502:mRNA:chr01.502");
		BlastResult br = p.getBlastResults().get("swissprot").get(0);
		assertEquals(2.0, br.getTokens().size(), 0.0);
		// test measurement of cumulative token-scores was triggered:
		assertTrue(p.getTokenScoreCalculator().getCumulativeTokenBitScores()
				.containsKey("dicer"));
		assertTrue(p.getTokenScoreCalculator()
				.getCumulativeTokenBlastDatabaseScores().containsKey("dicer"));
		assertTrue(p.getTokenScoreCalculator()
				.getCumulativeTokenOverlapScores().containsKey("dicer"));
		// test measurement of total token-scores was triggered:
		assertTrue(p.getTokenScoreCalculator().getTotalTokenBitScore() > 0.0);
		assertTrue(p.getTokenScoreCalculator()
				.getTotalTokenBlastDatabaseScore() > 0.0);
		assertTrue(p.getTokenScoreCalculator().getTotalTokenOverlapScore() > 0.0);
	}

	@Test
	public void testParseInterproResults() throws Exception {
		ahrd.setProteins(TestUtils.mockProteinDb());
		// TODO: The Interpro-Database should be mocked:
		InterproResult.initialiseInterproDb();
		ahrd.parseInterproResult();
		for (String iterProtAcc : ahrd.getProteins().keySet()) {
			Protein protein = ahrd.getProteins().get(iterProtAcc);
			assertTrue(protein.getInterproResults().size() > 0);
		}
	}

	@Test
	public void testParseGeneOntologyResults() throws Exception {
		ahrd.setProteins(TestUtils.mockProteinDb());
		ahrd.parseGeneOntologyResult();
		assertEquals(1, ahrd.getProteins().get("gene:chr01.502:mRNA:chr01.502")
				.getGoResults().size());
		assertEquals(2,
				ahrd.getProteins().get("gene:chr01.1056:mRNA:chr01.1056")
						.getGoResults().size());
	}

	//@Test
	public void testLoadBlastResultDomainAnnotationFromUniprotKB()
			throws InterruptedException {
		Protein p = TestUtils.mockProtein();
		List<BlastResult> brs = new ArrayList<BlastResult>();
		for (String accession : new String[] { "B5YXA4", "Q8XBZ3", "A7ZTQ8",
				"A8A6G3", "B1LL27", "B6I3T5", "B7L846", "B7M553", "B7MGC3",
				"B7N2E7", "B7NF19", "B7NR04", "B7UMG9", "C4ZYX9", "P03004",
				"Q1R4N5", "Q3YWB2", "B1IYP2", "B4F0U5", "P22837", "Q329B6",
				"B2TUS6", "Q7MQJ7", "Q8DDI9", "Q87TQ7", "Q31UV5", "A7N1E9",
				"P49996", "B5QUQ0", "Q8Z2N6", "A9MX77", "B4SYA8", "B4TAU9",
				"B4TN07", "B5BIL4", "B5EYX2", "P35891", "Q5PKU6", "C6DGH7",
				"B5FN10", "B5RFY7", "Q9KVX6", "P29440", "B7VGI4", "Q6CYR4",
				"A4TGL5", "A7FPB7", "A9R5R5", "Q1C0B8", "Q1CCJ8", "Q663T2",
				"Q8Z9U7", "B5XT51", "Q57I00", "Q7NAD3", "A3CYH5", "A6WH85",
				"A9KU72", "B8E3P4", "Q2NX49", "Q12TC8", "A1JT80", "C5BHC5",
				"A8G7Q2", "B1KCX3", "Q0HPD4", "Q8EKT2", "A1RDX7", "A4Y1A4",
				"A8FP46", "Q08A51", "Q0I0U8", "A0KR35", "B2VCE3", "A8GYE3",
				"Q5E8Z2", "B6EP46", "A3Q8S6", "B0TLA4", "B8CH71", "A1S1G9",
				"Q6LW50", "A0KEC3", "A4SH46", "C1DFU2", "C4L755", "Q4KKT0",
				"Q48AS7", "Q3KKG1", "B0B0A5", "Q500U7", "Q48QK0", "Q88BK3",
				"B0KEU9", "A5VWB8", "B1J3Y2", "P0A116", "P0A117", "Q1I2G4",
				"Q21PW4", "A5I9F1", "Q5X0L8", "Q5X990", "Q5ZZK8", "Q1R1P2",
				"B7V0N6", "Q02V80", "Q9I7C5", "A9KEU8", "A9N900", "B6J287",
				"B6J8S3", "Q83FD8", "A1TWJ0", "C5BKL9", "Q5L3Z2", "Q2SQZ9",
				"Q0ACS7", "A4IJ84", "A6VR65", "A9B496", "C5D327", "B7GFK8",
				"Q0VT30", "Q5QY39", "B2A2Y6", "Q65PM2", "B7HPR7", "B9IYG8",
				"Q73FK5", "B7JJB7", "C3LIC2", "C3P8P5", "Q63HG7", "Q6HQ03",
				"Q81W35", "C1ES08", "Q0TV64", "Q8XPG2", "A7GJR9", "B7HIH4",
				"Q81JD5", "B2THB4", "B7IS20", "P05648", "Q0SWX6", "B2UX43",
				"A7Z0C3", "A8F8Y4", "Q8EU88", "Q9RCA2", "P29434", "A1WWE0",
				"A6LPB1", "A0Q3U6", "Q7P259", "A9VM90", "Q31JS5", "Q1GXK1",
				"A5HXP7", "A7FPR6", "A5N457", "B9DXS7", "A7G9B0", "B1IDU3",
				"C1FPH3", "Q97N35", "C3KXQ7", "B1L1K6", "A0AEI7", "Q92FV2",
				"C1L2Z8", "B8DAQ9", "Q725H0", "Q8YAW2", "A1K1B4", "Q5WM31",
				"B0TAK8", "A3DHZ4", "Q5P4P0", "A8MEA0", "Q8RDL6", "Q4A180",
				"Q4LAL5", "Q2YD61", "Q8CQK7", "Q5HJZ9", "Q602N0", "C4KZZ3",
				"B1YGB2", "A4G154", "A6TJ76", "Q7VSE0", "Q7W2K5", "Q6GKU4",
				"B0K0W8", "B0KAG0", "Q82Y84", "A5INP2", "A6QD41", "A6TXF1",
				"A7WWM4", "A8YYS4", "P68865", "P68866", "P68867", "P68868",
				"Q2FKQ5", "Q2G2H5", "Q5HJZ5", "Q6GD89", "B8I3R2", "Q2YUP1",
				"A6STW2", "B8D6T1", "B8D8H7", "P57128", "C6E7Q5", "Q65VB8",
				"Q0I0Y7", "C1D6I2", "A4J0F0", "Q0AK27", "B5E7P6", "Q0KFR8",
				"A3N473", "Q3JXI6", "Q63YW5", "Q7WDJ9", "Q2KTI9", "B9E8Z7",
				"Q2STL6", "A1V7D9", "A3NPW7", "Q1LSI9", "A2S8D2", "A3MH48",
				"B2AFZ7", "B9DPX4", "B2JJ97", "Q67TK7", "A5TY69", "P49993",
				"P49991", "C1AIZ8", "Q39L82", "Q9CLQ4", "A9AI97", "B4E7D1",
				"A0K2M8", "B1YPZ0", "Q0BJW1", "A9HVA1", "B1K0Y8", "A7NFB8",
				"Q8XTV4", "Q47K71", "A5UP91", "B2UCF1", "B2SZ75", "Q3JF39",
				"A4J9S6", "Q147F0", "A6VK86", "P43742", "B2HI46", "A1U8S0",
				"A3PSD7", "Q1BG61", "A1T102", "A0PKB2", "Q9L7L7", "B8GBK7",
				"A0R7K1", "P0C557", "Q7VMW1", "A9WAN1", "B9LFG0", "P46388",
				"A5GDX1", "B1MDH6", "A1T0X4", "Q47U23", "Q9RYE7", "B9M7S1",
				"Q39ZS3", "A2SBM4", "P0A3A2", "P0A3A3", "B5YGT9", "B0CDM2",
				"Q890K8", "P49990", "A2BVM7", "A6W3V4", "Q74GG6", "Q839Z5",
				"Q03UE4", "Q6ABL5", "Q72H87", "Q9X9D5", "Q2JQW9", "B2J2B2",
				"C1CXJ1", "Q0I8T6", "A2BQ46", "Q8YVG9", "B9L0U6", "C3PE72",
				"Q9ZH75", "A2C122", "Q9ZH76", "B1Y547", "Q3MHA9", "B8HRT5",
				"Q8DL93", "P27902", "A0L3I7", "B1VPF0", "A8G3T0", "Q6AHN6",
				"Q31BW9", "Q7V6H4", "A3PBT9", "Q82FD8", "A5VHF3", "B2G4Y5",
				"A2C7X1", "Q3BZT1", "A5CLT3", "Q2JHS1", "Q8PRG2", "C0ZLE1",
				"A5GJL3", "B1X0X6", "B0TWF7", "Q5H715", "B2SUW3", "Q2P9M1",
				"C1B7S7", "B0RLI8", "Q8PEH5", "A9KPP1", "A4IVR6", "Q0BPB6",
				"A7N926", "Q5NIQ8", "B0RH69", "A0Q3U7", "B7K7Y7", "B2FUW1",
				"Q03I60", "Q38ZS4", "Q2A640", "Q5Z3Z8", "C5C7X4", "P21173",
				"B4SU11", "Q8NUD8", "Q0SAG7", "A4Q9R9", "A1AJX2", "B3W6N4",
				"B7JYF1", "Q51896", "Q8KGG6", "Q03D55", "Q6NKL7", "P49995",
				"Q5FN15", "Q8FUL7", "Q87FC6", "Q9PHE3", "Q04CX5", "Q1GC43",
				"Q4V0S8", "B0VAF3", "B2HZA7", "B7GUX5", "B7IBH7", "Q11AE3",
				"B0VMK0", "A3M0Q4", "C4K1R9", "A8YW41", "B3DP22", "Q8G6K0",
				"C3PP41", "Q1GKT2", "O87546", "A8GSY1", "B0BYF7", "Q92H56",
				"Q6FG21", "B7GSF9", "Q74M34", "Q4UMJ3", "A8EY77", "A8F284",
				"B0JGA6", "Q4JYF7", "Q7U605", "Q3AYH5", "A3PFL5", "Q3IY61",
				"Q68WD8", "Q59758", "Q5LWV4", "A8GP75", "Q1RI85", "A8GVN1",
				"B1XKQ0", "B3DWG6", "Q7NKK4", "Q28WI0", "Q5FUU1", "A9BEI9",
				"Q16DK6", "B2GEU8", "Q5FHH8", "Q5HBP0", "Q2GG27", "Q9S493",
				"Q5FAJ2", "A1KS02", "Q9JW45", "Q9JXS7", "Q07VS2", "A4YJ98",
				"A5E812", "Q13F98", "Q89W63", "B3QJZ9", "Q6NDV3", "Q3SMV0",
				"Q9CJJ2", "Q1QS94", "Q5PB48", "Q2J497", "A2RH74", "Q033I4",
				"P59567", "A6GYW8", "Q21DF6", "B8DV06", "A1QZM2", "Q5GT09",
				"Q73IZ0", "B2S0D9", "Q0SN72", "Q1J960", "C0MC62", "Q661I4",
				"C0M7C0", "B4U5G8", "B5E568", "C1C982", "O08397", "B1I6X3",
				"B2IQY9", "B8ZJH9", "C1CH81", "Q04N63", "Q8DRQ4", "C1CN66",
				"A2RBX3", "P0A3A4", "P0A3A6", "P0DA68", "P0DA69", "Q1JEA5",
				"Q1JJA7", "Q1JP60", "Q48VX2", "Q5XEM7", "Q8DWN9", "C1CU44",
				"B5RRN9", "B5RLZ3", "B9DSN7", "Q03N26", "Q5M226", "Q5M6L8",
				"Q3K425", "Q8E2I7", "Q8E7Z4", "P33768", "B7J203", "B5XIP2",
				"Q6MRS1", "P35890", "B9L735", "Q7MSY2", "Q058F9", "Q5L9D0",
				"Q64PL4", "Q823P0", "Q1MMD6", "Q83N52", "Q83NZ5", "Q8A5U5",
				"Q04HR6", "Q8G3E7", "Q57G10", "Q8YED5", "Q8UIH1", "A8EQT0",
				"B2S1V1", "O83047", "Q6G526", "Q9Z8B9", "Q98BG9", "Q6G0V6",
				"A0RLX8", "O84277", "Q9PKB9", "A7GVR3", "Q30UP9", "A7IB67",
				"B2KAM4", "Q6MEG8", "A7ZAW7", "A9NE65", "B2RGM5", "Q7MXZ1",
				"A8FJG5", "Q9PJB0", "Q9KHU8", "A7H194", "Q5HXF5", "A1VX79",
				"B8GWW5", "P0CAU4", "A5IIK7", "P46798", "B1LC08", "Q7VH57",
				"A8F346", "B0T135", "Q1GN68", "Q6YRL3", "C5CH91", "B3R0L6",
				"Q2NKC5", "Q30YX9", "B8FFL8", "Q1MSG8", "B8DN19", "Q6ARL8",
				"Q729U6", "Q17ZQ6", "A0LE53", "B5Z9E3", "P34028", "B2UVS4",
				"Q1CRG8", "Q9ZJ96", "B6JP23", "O26057", "Q72WD6", "Q8FA34",
				"Q04WF7", "Q056V2", "P24116", "A6LIY2", "B0S907", "B0SK31",
				"Q6MUM7", "B7IF65", "Q823E6", "Q9Z8M9", "Q6MC93", "Q9PKE4",
				"O84252", "O66659", "Q6F2A9", "Q59549", "P35888", "B1AHY7",
				"Q9PRE2", "P35892", "P35889", "P35907", "B4EY85", "A6TCA8",
				"B5XNQ6", "B7LKE6", "A1ADY9", "A7ZPT9", "A8A2Y8", "B1IWG4",
				"B1LNE6", "B1XAX1", "B2TXS1", "B5Z033", "B6I568", "B7LCN6",
				"B7M7K0", "B7MHX7", "B7MYC3", "B7N678", "B7NQN5", "B7UGN4",
				"C4ZX71", "P69931", "P69932", "P69933", "P69934", "Q0T224",
				"Q0TEZ2", "Q1R8P2", "Q31XZ6", "Q32D72", "Q3YZ56", "A1JL00",
				"A4TMP0", "B1JSF7", "B2K9J7", "Q1C5P5", "Q1CK37", "Q668E8",
				"A7FG40", "A9QZX0", "Q8ZCC1", "A9MHP3", "A9N2Y3", "B4T0M5",
				"B4TD71", "B4TR72", "B5F171", "B5FQI8", "B5R558", "B5RCW8",
				"C0PYR0", "Q57LL3", "Q7CQ21", "Q8XEQ0", "B5BB08", "Q5PL41",
				"O86235", "A8AD95", "Q6D7R8", "Q7N3G5", "A7ML17", "P77546" }) {
			brs.add(new BlastResult(accession, 0.01, "description", 10, 20, 10,
					20, 200, 10.0, "swissprot"));
		}
		p.getBlastResults().put("swissprot", brs);

		// Test:
		this.ahrd.loadBlastResultDomainAnnotationFromUniprotKB(p);
		assertEquals(8, DomainScoreCalculator
				.getBlastResultAccessionsToInterproIds().get("B5YXA4").size());
	}

}
