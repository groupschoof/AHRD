package ahrd.view;

import static ahrd.controller.Settings.getSettings;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import ahrd.model.Protein;

import com.google.gson.Gson;

public class JsonOutputWriter extends AbstractOutputWriter {

	public List<ProteinOutput> proteinOutputs = new ArrayList<ProteinOutput>();

	public JsonOutputWriter(Collection<Protein> proteins) {
		super(proteins);
		for (Protein prot : getProteins()) {
			this.proteinOutputs.add(new ProteinOutput(prot));
		}
	}

	@Override
	public void writeOutput() throws IOException {
		BufferedWriter o = new BufferedWriter(new FileWriter(getSettings()
				.getPathToOutput()));
		Gson gson = new Gson();
		o.write(gson.toJson(this.proteinOutputs) + "\n");
		o.close();
	}
}
