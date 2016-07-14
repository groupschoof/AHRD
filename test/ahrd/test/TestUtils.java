package ahrd.test;

import static ahrd.controller.Settings.setSettings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import ahrd.controller.Settings;
import ahrd.model.BlastResult;
import ahrd.model.LexicalScoreCalculator;
import ahrd.model.Protein;

public class TestUtils {

	// Mock Classes:
	protected static class LexicalScoreCalculatorMock extends LexicalScoreCalculator {

		public LexicalScoreCalculatorMock(Protein protein) {
			super(protein);
		}

		public double lexicalScore(BlastResult br) {
			return 0.70;
		}
	}

	protected static class LexicalScoreCalculatorFixedGoScoreMock extends LexicalScoreCalculator {

		public LexicalScoreCalculatorFixedGoScoreMock(Protein protein) {
			super(protein);
		}

		public double geneOntologyScore(BlastResult blastResult) {
			return 1.11;
		}
	}

	// END Mock Classes!

	/**
	 * For each Test-Class (Suite) the AHRD-Settings are set by calling this
	 * method with the slight disadvantage of the repeated reading out the test
	 * ahrd_input.yml.
	 * 
	 * @throws IOException
	 */
	public static void initTestSettings() throws IOException {
		setSettings(new Settings("./test/resources/ahrd_input.yml"));
	}

	public static BlastResult mockBlastResult(String acc, Double eValue, String descLine, int queryStart, int queryEnd,
			int subjectStart, int subjectEnd, int subjectLength, Double bitScore, String dbName, Set<String> tokens) {
		BlastResult br = new BlastResult(acc, eValue, descLine, queryStart, queryEnd, subjectStart, subjectEnd,
				subjectLength, bitScore, dbName);
		br.setShortAccession(acc);
		br.setTokens(tokens);
		return br;
	}

	public static Map<String, Protein> mockProteinDb() {
		// TODO: Actually mock the Protein-DB!
		Map<String, Protein> protDb = new HashMap<String, Protein>();
		Protein prot = new Protein("gene:chr01.1056:mRNA:chr01.1056", 892);
		protDb.put(prot.getAccession(), prot);
		prot = new Protein("gene:chr01.502:mRNA:chr01.502", 108);
		protDb.put(prot.getAccession(), prot);
		return protDb;
	}

	public static List<BlastResult> mockBlastResults() {
		List<BlastResult> blastResults = new ArrayList<BlastResult>();
		blastResults.add(new BlastResult("accession_1", 1.0, "description One", 10, 20, 10, 20, 200, 30, "swissprot"));
		blastResults.add(new BlastResult("accession_2", 2.0, "description Two", 10, 20, 10, 20, 200, 30, "swissprot"));
		blastResults.add(new BlastResult("accession_3", 3.0, "Putative - sUbFaMilY;, \" activity|, bad", 10, 20, 10, 20,
				200, 30, "swissprot"));
		blastResults.add(new BlastResult("accession_4", 4.0, "family subfamily activity NADH-Dehydrogenase", 10, 20, 10,
				20, 200, 30, "swissprot"));
		blastResults.add(new BlastResult("accession_5", 5.0, "description_5 Fly-Wing formation", 10, 20, 10, 20, 200,
				30, "swissprot"));
		return blastResults;
	}

	public static List<BlastResult> mockBlastResultsForDescCalcTest() {
		List<BlastResult> blastResults = new ArrayList<BlastResult>();
		blastResults.add(mockBlastResult("accession_1", 1.0, "description One", 10, 20, 10, 20, 200, 30.0, "swissprot",
				new HashSet<String>(Arrays.asList("description", "one"))));
		blastResults.add(mockBlastResult("accession_2", 2.0, "description Two", 10, 20, 10, 20, 200, 30.0, "swissprot",
				new HashSet<String>(Arrays.asList("description", "two"))));
		blastResults.add(mockBlastResult("accession_3", 3.0, "Putative - sUbFaMilY;, \" activity|, bad", 10, 20, 10, 20,
				200, 30.0, "swissprot",
				new HashSet<String>(Arrays.asList("putative", "subfamily", "activity", "bad"))));
		blastResults.add(mockBlastResult("accession_4", 4.0, "family subfamily activity NADH-Dehydrogenase", 10, 20, 10,
				20, 200, 30.0, "swissprot",
				new HashSet<String>(Arrays.asList("family", "subfamily", "activity", "nadh", "dehydrogenase"))));
		return blastResults;
	}

	public static Map<String, Double> mockTokenScoresForDescCalcTest() {
		Map<String, Double> tknScrs = new HashMap<String, Double>();
		tknScrs.put("description", 0.5);
		tknScrs.put("One", 1.0);
		tknScrs.put("Two", 1.5);
		tknScrs.put("Putative", 2.0);
		tknScrs.put("sUbFaMilY", 2.5);
		tknScrs.put("activity", 3.0);
		tknScrs.put("bad", 3.25);
		tknScrs.put("family", 3.3);
		tknScrs.put("subfamily", 3.4);
		tknScrs.put("NADH", 3.5);
		tknScrs.put("Dehydrogenase", 4.0);
		tknScrs.put("5", 4.5);
		tknScrs.put("Fly", 5.0);
		tknScrs.put("Wing", 5.5);
		tknScrs.put("formation", 6.0);
		return tknScrs;
	}

	public static Map<String, Integer> mockDescriptionLineFrequenciesForDescCalcTest() {
		Map<String, Integer> descLineFreaks = new HashMap<String, Integer>();
		descLineFreaks.put("descriptionone", 1);
		descLineFreaks.put("descriptiontwo", 2);
		descLineFreaks.put("activitybadputativesubfamily", 3);
		descLineFreaks.put("activitydehydrogenasefamilynadhsubfamily", 4);
		descLineFreaks.put("5descriptionflyformationwing", 5);
		return descLineFreaks;
	}

	public static List<BlastResult> mockBlastResultsWithTokens() {
		List<BlastResult> blastResults = new ArrayList<BlastResult>();
		BlastResult one = new BlastResult("accession_1", 1.0, "one two", 10, 20, 10, 20, 200, 30, "swissprot");
		BlastResult two = new BlastResult("accession_2", 2.0, "three", 10, 20, 10, 20, 200, 30, "swissprot");
		String elements1[] = { "one", "two" };
		String elements2[] = { "three" };
		Set<String> tokens1 = new HashSet<String>(Arrays.asList(elements1));
		Set<String> tokens2 = new HashSet<String>(Arrays.asList(elements2));

		one.setTokens(new HashSet<String>(tokens1));
		two.setTokens(new HashSet<String>(tokens2));
		blastResults.add(one);
		blastResults.add(two);

		return blastResults;
	}

	public static void mockCumulativeTokenScores(Protein p, String token, double mockBase) {
		p.getTokenScoreCalculator().getCumulativeTokenBitScores().put(token, 5 * mockBase);
		p.getTokenScoreCalculator().getCumulativeTokenBlastDatabaseScores().put(token, 10 * mockBase);
		p.getTokenScoreCalculator().getCumulativeTokenOverlapScores().put(token, 0.05 * mockBase);
	}

	public static BlastResult mockBlastResult() {
		BlastResult br = new BlastResult("accession_1", 1.0, "one two three", 10, 20, 10, 20, 200, 30, "swissprot");
		br.getTokens().add("one");
		br.getTokens().add("two");
		br.getTokens().add("three");
		return br;
	}

	/**
	 * Mocks Protein of length 200.
	 * 
	 * @return new Protein
	 */
	public static Protein mockProtein() {
		return new Protein("sweet_sheep_protein", 200);
	}

	/**
	 * Used in test functions in class DescriptionScoreCalculatorTest
	 * 
	 * @return new Protein
	 */
	public static Protein mockProteinAndBlastResultsForDescriptionScoreCalculatorTest() {
		Protein p = mockProtein();
		// Sprot
		p.getBlastResults().put("swissprot", TestUtils.mockBlastResultsForDescCalcTest());
		// trEMBL
		p.getBlastResults().put("trembl",
				Arrays.asList(mockBlastResult("accession_5", 5.0, "description_5 Fly-Wing formation", 10, 20, 10, 20,
						200, 30.0, "trembl",
						new HashSet<String>(Arrays.asList("description", "5", "fly", "wing", "formation")))));
		p.setLexicalScoreCalculator(new LexicalScoreCalculatorMock(p));
		p.getDescriptionScoreCalculator().setMaxBitScore(30.0);
		return p;
	}

	public static Map<String, Set<String>> mockReferenceGoAnnotationsForDescriptionScoreCalculatorTest() {
		Map<String, Set<String>> refGos = new HashMap<String, Set<String>>();
		refGos.put("accession_1", new HashSet<String>(Arrays.asList("GO:1234567", "GO:7654321")));
		refGos.put("accession_4", new HashSet<String>(Arrays.asList("GO:1726354", "GO:7162534")));
		return refGos;
	}
}
