package ahrd.test;

import static ahrd.controller.AddEuclideanDistsToPathLogTable.isAccepted;
import static ahrd.controller.AddEuclideanDistsToPathLogTable.main;
import static ahrd.controller.AddEuclideanDistsToPathLogTable.measureEuclideanDist;
import static ahrd.controller.AddEuclideanDistsToPathLogTable.normalizeDatabaseWeights;
import static ahrd.controller.AddEuclideanDistsToPathLogTable.parseCurrentPathLogLine;
import static ahrd.controller.AddEuclideanDistsToPathLogTable.parseDatabaseWeightColumnArg;
import static ahrd.controller.Utils.fromFile;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;

import org.junit.Test;

public class AddEuclideanDistsToPathLogTableTest {

	@Test
	public void testMeasureEuclideanDist() {
		double[] x = new double[] { 0, 0, 0, 0 };
		double[] y = new double[] { 1, 1, 1, 1 };
		assertEquals(2.0, measureEuclideanDist(x, y), 0.0);
		x = new double[] { 0.10000, 0.79625, 0.11375, 0.09000, 0.30000,
				30.00000, 0.10000, 30.00000 };
		y = new double[] { 0.10000, 0.79625, 0.11375, 0.09000, 0.30000,
				30.00000, 0.10000, 30.27175 };
		assertEquals(0.27175, measureEuclideanDist(x, y), 0.00001);
	}

	@Test
	public void testIsAccepted() {
		assertTrue(isAccepted(
				"50000\t0.5263841932366594\t0.0\t3\t0.524761\t0.157966\t0.1\t0.79625\t0.11375\t0.09\t30\t30\t10\t30",
				3));
		assertTrue(!isAccepted(
				"49999\t0.5223735235001833\t-0.004010669736476125\t0\t0.518299\t0.157462\t0.1\t0.79625\t0.11375\t0.09\t30\t30\t10\t30.271748",
				3));
	}

	@Test
	public void testNormalizeDatabaseWeights() {
		double[] x = new double[] { 0.10000, 0.79625, 0.11375, 0.09000,
				0.30000, 30.00000, 0.10000, 30.00000 };
		double[] expctd = new double[] { 0.10000, 0.79625, 0.11375, 0.09000,
				0.30000, 0.3, 0.10000, 0.3 };
		double[] actual = normalizeDatabaseWeights(x, new int[] { 5, 7 });
		for (int i = 0; i < expctd.length; i++) {
			assertEquals(expctd[i], actual[i], 0.000001);
		}
	}

	@Test
	public void testParseCurrentPathLogLine() {
		double[] p = parseCurrentPathLogLine(
				"49999\t0.5223735235001833\t-0.004010669736476125\t0\t0.518299\t0.157462\t0.1\t0.79625\t0.11375\t0.09\t30\t30\t10\t30.2710748",
				6, 14);
		double[] expct = new double[] { 0.1, 0.79625, 0.11375, 0.09, 30, 30,
				10, 30.2710748 };
		assertEquals(expct.length, p.length);
		for (int i = 0; i < expct.length; i++) {
			assertEquals(expct[i], p[i], 0.000001);
		}
	}

	@Test
	public void testParseDatabaseWeightColumnArg() {
		int[] dbWghtCols = parseDatabaseWeightColumnArg("2,4,6", 1);
		assertEquals(3, dbWghtCols.length);
		assertEquals(1, dbWghtCols[0]);
		assertEquals(3, dbWghtCols[1]);
		assertEquals(5, dbWghtCols[2]);
	}

	@Test
	public void testMain() throws IOException {
		try {
			main(new String[] { "test/resources/path_log_tbl.tsv", "6", "13",
					"3", "10,12", "./test/resources/tmp_path_log.tsv" });
			List<String> lines = fromFile("./test/resources/tmp_path_log.tsv");
			assertEquals(30, lines.size());
			// Check the first lines:
			assertEquals(
					"Temperature\tAverage Evaluation-Score(F-Score)\tDiff-to-curr-Accepted\tAccepted\tAverage True-Positive-Rate\tAverage False-Positive-Rate\tDescription-Score-Relative-Description-Frequency-Weight\tToken-Score-Bit-Score-Weight\tToken-Score-Database-Score-Weight\tToken-Score-Overlap-Score-Weight\tswissprot-Weight\tswissprot-Description-Score-Bit-Score-Weight\ttrembl-Weight\ttrembl-Description-Score-Bit-Score-Weight\tEuclidean.Distance.2.Curr.Acptd",
					lines.get(0));
			assertEquals(
					"50000\t0.5263841932366594\t0.0\t3\t0.524761\t0.157966\t0.1\t0.79625\t0.11375\t0.09\t30\t30\t10\t30\t0.0",
					lines.get(1));
			assertEquals(
					"49999\t0.5223735235001833\t-0.004010669736476125\t1\t0.518299\t0.157462\t0.1\t0.79625\t0.11375\t0.09\t30\t30\t10\t30.271748\t0.27174799999999877",
					lines.get(2));
			assertEquals(
					"49998\t0.5251306827089552\t-0.0012535105277041714\t0\t0.522247\t0.156624\t0.1\t0.79625\t0.11375\t0.09\t30\t30\t10\t29.491517\t0.780230999999997",
					lines.get(3));
		} finally {
			// clean up:
			Files.deleteIfExists(Paths.get("./test/resources/tmp_path_log.tsv"));
		}
	}
}
