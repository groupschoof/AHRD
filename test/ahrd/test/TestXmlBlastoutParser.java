package ahrd.test;

import java.io.File;
import java.util.Collections;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.Unmarshaller;

import org.junit.Test;

import ahrd.model.blast.BlastOutput;
import ahrd.model.blast.Hit;
import ahrd.model.blast.Hsp;
import ahrd.model.blast.Iteration;
import ahrd.model.blast.IterationHits;

public class TestXmlBlastoutParser {

	/**
	 * If you want to test quickly from the command line with any custom Blast
	 * XML output file
	 * 
	 * @param args
	 */
	public static void main(String[] args) {
		try {
			parse(args[0]);
		} catch (Exception e) {
			e.printStackTrace(System.err);
		}
	}

	@Test
	public void test() throws Exception {
		parse("./test/resources/blast_output.xml");
	}

	public static void parse(String path) throws Exception {
		JAXBContext jc = JAXBContext.newInstance(BlastOutput.class);
		Unmarshaller u = jc.createUnmarshaller();
		BlastOutput bo = (BlastOutput) u.unmarshal(new File(path));

		for (Iteration iter : bo.getBlastOutputIterations().getIteration()) {
			System.out.println("QueryDef: " + iter.getIterationQueryDef());
			System.out.println("QueryId: " + iter.getIterationQueryID());
			IterationHits iterHits = iter.getIterationHits();
			if (iterHits != null) {
				for (Hit hit : iterHits.getHit()) {
					System.out.println("Hit-Description: " + hit.getHitDef());
					System.out.println("Hit-Accession: "
							+ hit.getHitAccession());
					Hsp bestHsp = Collections.min(hit.getHitHsps().getHsp());
					if (bestHsp != null) {
						System.out.println("Best HSP E-Value: "
								+ bestHsp.getHspEvalue());
						System.out.println("Best HSP Hit-From: "
								+ bestHsp.getHspHitFrom());
						System.out.println("Best HSP Hit-To: "
								+ bestHsp.getHspHitTo());
						System.out.println("Best HSP Query-From: "
								+ bestHsp.getHspQueryFrom());
						System.out.println("Best HSP Query-To: "
								+ bestHsp.getHspQueryTo());
					} else {
						System.out.println("No best HSP found");
					}
				}
			} else {
				System.out.println(">>>>>> NO HITS <<<<<<<<<");
			}

		}
		return;
	}
}
