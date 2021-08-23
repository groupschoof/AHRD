package ahrd.test;

import static ahrd.controller.Settings.getSettings;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.junit.Before;
import org.junit.Test;

import ahrd.controller.DescriptionParameters;
import ahrd.controller.Parameters;
import ahrd.controller.Settings;

public class DescriptionParametersTest {

	@Before
	public void setUp() throws IOException {
		TestUtils.initTestSettings();
	}

	@Test
	public void testParameterToMutateRandomIndex() {
		Parameters p = getSettings().getDescriptionParameters();
		Set<Integer> inds = new HashSet<>();
		for (int i = 0; i < 500; i++) {
			inds.add(p.parameterToMutateRandomIndex());
		}
		assertEquals(p.numberOfNonDbParameters() + 2 * 3, inds.size());
		// Parameter-Indices: 0..2 + 2 * 3 (#Blast-Databases) = 8
		// Indices 0 to 8 should all be present:
		for (int r = 0; r <= 8; r++) {
			assertTrue(
					"Parameter-Index "
							+ r
							+ " should be present in Set of randomly selected Parameter-Indices.",
					inds.contains(r));
		}
	}

	@Test
	public void testNormalizeTokenScoreWeights() {
		Parameters p = getSettings().getDescriptionParameters();
		p.setTokenScoreBitScoreWeight(0.5);
		p.setTokenScoreDatabaseScoreWeight(0.75);
		p.setTokenScoreOverlapScoreWeight(1.0);
		p.normalizeTokenScoreWeights();
		assertEquals(
				"The three weights in the Token-Score-Formula should sum up to 1.0",
				1.0,
				p.getTokenScoreBitScoreWeight()
						+ p.getTokenScoreDatabaseScoreWeight()
						+ p.getTokenScoreOverlapScoreWeight(), 0.0001);
		List<String> sortedDistinctBlastDatabaseNames = new ArrayList<String>();
		sortedDistinctBlastDatabaseNames.addAll(getSettings().getBlastDatabases());
		Collections.sort(sortedDistinctBlastDatabaseNames);
		p = p.randomParameters(sortedDistinctBlastDatabaseNames);
		assertEquals(
				"The three randomly generated and normalized weights in the Token-Score-Formula should sum up to 1.0",
				1.0,
				p.getTokenScoreBitScoreWeight()	+ p.getTokenScoreDatabaseScoreWeight() + p.getTokenScoreOverlapScoreWeight(), 0.001);
		p.mutateTokenScoreBitScoreWeight();
		assertEquals(
				"The mutated TokenScoreBitScoreWeight and the other token score weights in the Token-Score-Formula should sum up to 1.0 after beeing normalized",
				1.0,
				p.getTokenScoreBitScoreWeight()	+ p.getTokenScoreDatabaseScoreWeight() + p.getTokenScoreOverlapScoreWeight(), 0.001);
		p.mutateTokenScoreDatabaseScoreWeight();
		assertEquals(
				"The mutated TokenScoreDatabaseScoreWeight and the other token score weights in the Token-Score-Formula should sum up to 1.0 after beeing normalized",
				1.0,
				p.getTokenScoreBitScoreWeight()	+ p.getTokenScoreDatabaseScoreWeight() + p.getTokenScoreOverlapScoreWeight(), 0.001);
		p.mutateTokenScoreOverlapScoreWeight();
		assertEquals(
				"The mutated TokenScoreOverlapScoreWeight and the other token score weights in the Token-Score-Formula should sum up to 1.0 after beeing normalized",
				1.0,
				p.getTokenScoreBitScoreWeight()	+ p.getTokenScoreDatabaseScoreWeight() + p.getTokenScoreOverlapScoreWeight(), 0.001);
	}

	@Test
	public void testMutateBlastDbWeight() {
		Parameters p = getSettings().getDescriptionParameters();
		Map<String, Integer> origDbWeights = new HashMap<String, Integer>();
		for (String blastDB : p.getBlastDatabases()) {
			origDbWeights
					.put(blastDB, p.getBlastDbWeight(blastDB));
		}
		// test:
		for (String blastDB : p.getBlastDatabases()) {
			p.mutateBlastDatabaseWeight(blastDB);
			assertTrue(
					"Blast-Database-Weights of Blast-Database " + blastDB
							+ " should be mutated",
					!origDbWeights.get(blastDB).equals(
							p.getBlastDbWeight(blastDB)));
		}
	}

	@Test
	public void testRandomBlastDbName() {
		Parameters p = getSettings().getDescriptionParameters();
		List<String> randBlastDbs = new ArrayList<String>();
		for (int i = 0; i < 1000000; i++) {
			String r = p.randomBlastDatabaseName();
			assertTrue(
					"A random Blast-Database-Name must be one of the original ones.",
					p.getBlastDatabases().contains(r));
			randBlastDbs.add(r);
		}
		for (String blastDB : p.getBlastDatabases()) {
			assertTrue(
					"Randomly selecting a Blast-Database 1,000,000 times should contain Blast-Database '"
							+ blastDB + "'.", randBlastDbs.contains(blastDB));
		}
	}

	@Test
	public void testMutateAnnotationScoreBitScoreWeight() {
		Parameters p = getSettings().getDescriptionParameters();
		Map<String, Double> origDsBsWs = new HashMap<String, Double>();
		for (String blastDB : p.getBlastDatabases()) {
			origDsBsWs.put(blastDB, p.getAnnotationScoreBitScoreWeight(blastDB));
		}
		// test:
		for (String blastDB : p.getBlastDatabases()) {
			p.mutateAnnotationScoreBitScoreWeight(blastDB);
			assertTrue(
					"Blast-Database-Description-Score-Bit-Score-Weights of Blast-DB "
							+ blastDB + " should have been mutated.",
					!origDsBsWs.get(blastDB).equals(
							p.getAnnotationScoreBitScoreWeight(blastDB)));
		}
	}

	@Test
	public void testMutateTokenScoreBitScoreWeight() {
		Parameters p = getSettings().getDescriptionParameters();
		Double dbsw = p.getTokenScoreDatabaseScoreWeight();
		Double osw = p.getTokenScoreOverlapScoreWeight();
		Double bsw = p.getTokenScoreBitScoreWeight();
		// test:
		p.mutateTokenScoreBitScoreWeight();
		assertEquals(
				1.0,
				p.getTokenScoreBitScoreWeight()
						+ p.getTokenScoreDatabaseScoreWeight()
						+ p.getTokenScoreOverlapScoreWeight(), 0.001);
		assertTrue(
				"All three Token-Score-Weights should have changed.",
				!dbsw.equals(p.getTokenScoreDatabaseScoreWeight())
						&& !osw.equals(p.getTokenScoreOverlapScoreWeight())
						&& !bsw.equals(p.getTokenScoreBitScoreWeight()));
	}
	
	@Test
	public void testMutateTokenScoreOverlapScoreWeight() {
		Parameters p = getSettings().getDescriptionParameters();
		Double dbsw = p.getTokenScoreDatabaseScoreWeight();
		Double osw = p.getTokenScoreOverlapScoreWeight();
		Double bsw = p.getTokenScoreBitScoreWeight();
		// test:
		p.mutateTokenScoreOverlapScoreWeight();
		assertEquals(
				1.0,
				p.getTokenScoreBitScoreWeight()
						+ p.getTokenScoreDatabaseScoreWeight()
						+ p.getTokenScoreOverlapScoreWeight(), 0.001);
		assertTrue(
				"All three Token-Score-Weights should have changed.",
				!dbsw.equals(p.getTokenScoreDatabaseScoreWeight())
						&& !osw.equals(p.getTokenScoreOverlapScoreWeight())
						&& !bsw.equals(p.getTokenScoreBitScoreWeight()));
	}

	@Test
	public void testMutateTokenScoreDatabaseScoreWeight() {
		Parameters p = getSettings().getDescriptionParameters();
		Double dbsw = p.getTokenScoreDatabaseScoreWeight();
		Double osw = p.getTokenScoreOverlapScoreWeight();
		Double bsw = p.getTokenScoreBitScoreWeight();
		// test:
		p.mutateTokenScoreDatabaseScoreWeight();
		assertEquals(
				1.0,
				p.getTokenScoreBitScoreWeight()
						+ p.getTokenScoreDatabaseScoreWeight()
						+ p.getTokenScoreOverlapScoreWeight(), 0.001);
		assertTrue(
				"All three Token-Score-Weights should have changed.",
				!dbsw.equals(p.getTokenScoreDatabaseScoreWeight())
						&& !osw.equals(p.getTokenScoreOverlapScoreWeight())
						&& !bsw.equals(p.getTokenScoreBitScoreWeight()));
	}
	
	@Test
	public void testNeighbour() {
		// test
		DescriptionParameters n = getSettings().getDescriptionParameters().neighbour(0.0);
		assertNotNull("The neighbor of current Settings must not be NULL.", n);
		assertTrue(
				"The neighbour of current Settings must not be the same object as Settings. - Expecting a slightly changed CLONE!",
				!n.equals(getSettings().getDescriptionParameters()));
		// assert "slight difference" to currently set Parameters:
		Settings s = getSettings();
		boolean blastParamDiff = false;
		for (String blastDbName : getSettings().getBlastDatabases()) {
			blastParamDiff = blastParamDiff
					|| (!n.getBlastDbWeight(blastDbName).equals(s.getDescriptionBlastDbWeight(blastDbName)) 
					|| !n.getAnnotationScoreBitScoreWeight(blastDbName).equals(s.getDescriptionScoreBitScoreWeight(blastDbName)));
		}
		assertTrue(
				"The cloned neighbour must differ in exactly one AHRD-parameter-field from the currently set Settings.",
				(!n.getTokenScoreBitScoreWeight().equals(s.getDescriptionTokenScoreBitScoreWeight())
						|| !n.getTokenScoreDatabaseScoreWeight().equals(s.getDescriptionTokenScoreDatabaseScoreWeight()) 
						|| !n.getTokenScoreOverlapScoreWeight().equals(s.getDescriptionTokenScoreOverlapScoreWeight()))
						|| blastParamDiff);
		// Extreme Score-Increase should result in mutation of the same
		// parameter:
		for (int paramInd = 0; paramInd < n.numberOfNonDbParameters() + 2 * getSettings().getSortedBlastDatabases().size(); paramInd++) {
			n.setLastMutatedParameter(paramInd);
			DescriptionParameters n2 = n.neighbour(1.0);
			assertEquals(
					"The neighbor must remember which Parameter has been mutated to evolve him from its parent.",
					Integer.valueOf(paramInd), n2.getLastMutatedParameter());
			if (paramInd < n.numberOfNonDbParameters()) {
				switch (paramInd) {
				case 0: assertTrue("TokenScoreBitScoreWeight should have been mutated.",!n.getTokenScoreBitScoreWeight().equals(n2.getTokenScoreBitScoreWeight())); break;
				case 1: assertTrue("TokenScoreDatabaseScoreWeight should have been mutated.",!n.getTokenScoreDatabaseScoreWeight().equals(n2.getTokenScoreDatabaseScoreWeight())); break; 
				case 2: assertTrue("TokenScoreOverlapScoreWeight should have been mutated.",!n.getTokenScoreOverlapScoreWeight().equals(n2.getTokenScoreOverlapScoreWeight())); break;  
				}
			} else {
				String blastDbName = getSettings().getSortedBlastDatabases()
						.get((int) ((paramInd - n.numberOfNonDbParameters()) / 2.0));
				boolean mutatedBlastDbWeight = ! ((paramInd - n.numberOfNonDbParameters()) % 2 == 1);
				if (mutatedBlastDbWeight)
					assertTrue(
							"BlastDatabaseWeight of db " + blastDbName
									+ " should have been mutated.",
							!n.getBlastDbWeight(blastDbName).equals(
									n2.getBlastDbWeight(blastDbName)));
				else
					assertTrue(
							"DescriptionScoreBitScoreWeight of db "
									+ blastDbName
									+ " should have been mutated.",
							!n.getAnnotationScoreBitScoreWeight(blastDbName)
									.equals(n2
											.getAnnotationScoreBitScoreWeight(blastDbName)));
			}
		}
	}

	@Test
	public void testClone() {
		Parameters p = getSettings().getDescriptionParameters();
		p.setAvgEvaluationScore(0.8);
		Parameters c = p.clone();
		// test
		assertTrue("A clone should not be it's 'parent'.", p != c);
		// All primitive types get automatically cloned by Object.clone(),
		// so just test the Maps:
		for (String blastDb : p.getBlastDatabases()) {
			// For some very strange reason, calls to the getter
			// 'getBlastDbWeight(String blastDatabaseName)' return
			// Integer-Object, as in this case Integer.parseInt returns the same
			// object for the equal int-values. Maybe this is some optimization
			// of Memory-Usage in the JVM? - Anyway, this is the reason, why
			// here we compare the original Strings, held in the Map
			// blastDbSettings:
			assertTrue(
					"The BlastDatabaseWeights should not be the same objects.",
					(p.getBlastDbParameters().get(blastDb)
							.get(Settings.BLAST_DB_WEIGHT_KEY) != c
							.getBlastDbParameters().get(blastDb)
							.get(Settings.BLAST_DB_WEIGHT_KEY)));
			assertTrue(
					"The BlastDatabase's DescriptionScoreBitScoreWeight should not be the same objects.",
					System.identityHashCode(p
							.getAnnotationScoreBitScoreWeight(blastDb)) != System.identityHashCode(c
							.getAnnotationScoreBitScoreWeight(blastDb)));
		}
		// Test passing on the average evaluation score:
		assertEquals(p.getAvgEvaluationScore(), c.getAvgEvaluationScore(), 0.0);
		// Assure they are different Objects. As the operator != does not reveal
		// this, we set the clone to a different value than its parent:
		c.setAvgEvaluationScore(0.7);
		assertTrue(
				"Cloning should result in the average evaluation scores being different Objects.",
				p.getAvgEvaluationScore() != c.getAvgEvaluationScore());
	}

	@Test
	public void testEqualityAndHashCode() {
		Parameters p = getSettings().getDescriptionParameters();
		// Also test evaluation-scores:
		p.setAvgEvaluationScore(100.00);
		p.setAvgPrecision(0.5);
		p.setAvgRecall(0.5);
		Parameters c = p.clone();
		assertTrue("Cloned Parameters should still be equal.", p.equals(c));
		Set<Parameters> h = new HashSet<>();
		h.add(p);
		assertTrue(
				"Equal Parameters should not be added to a mathematical Set twice.",
				!h.add(c));
		// Changing scores should still fulfill equality:
		c.setAvgEvaluationScore(123.45);
		assertTrue(!p.getAvgEvaluationScore().equals(c.getAvgEvaluationScore()));
		assertTrue("Changing scores should still fulfill equality.",
				p.equals(c));
		// Change clone and check, if it is now unequal to parent:
		c.setTokenScoreBitScoreWeight(123.00);
		assertTrue("Changed clone should not be equal to parent", !p.equals(c));
		assertTrue("Unequal Parameters should be added to a mathematical Set.",
				h.add(c));
	}

	@Test
	public void testEquals() {
		Parameters p = getSettings().getDescriptionParameters();
		Parameters q = p.clone();
		assertEquals("Clones should be equal", p, q);
		q.setTokenScoreBitScoreWeight(1.5);
		assertTrue(
				"Changed TokenScoreBitScoreWeight should result in inequality.",
				!p.equals(q));
		q = p.clone();
		q.setTokenScoreDatabaseScoreWeight(1.5);
		assertTrue(
				"Changed TokenScoreDatabaseScoreWeight should result in inequality.",
				!p.equals(q));
		q = p.clone();
		q.setTokenScoreOverlapScoreWeight(1.5);
		assertTrue(
				"Changed TokenScoreOverlapScoreWeight should result in inequality.",
				!p.equals(q));
		q = p.clone();
		q.setBlastDbWeight("swissprot", "696");
		assertTrue(
				"Changed Swissprot-BlastDbWeight should result in inequality.",
				!p.equals(q));
		q = p.clone();
		q.setAnnotationScoreBitScoreWeight("swissprot", "1.5");
		assertTrue(
				"Changed Swissprot-AnnotationScoreBitScoreWeight should result in inequality.",
				!p.equals(q));
	}

}
