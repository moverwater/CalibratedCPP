package calibratedcpp.beauti;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.List;

import beast.base.core.BEASTInterface;
import beast.base.core.Input;
import beastfx.app.beauti.TreeDistributionInputEditor;
import beastfx.app.inputeditor.BeautiDoc;
import beastfx.app.util.FXUtils;
import calibratedcpp.CalibratedCoalescentPointProcess;
import javafx.scene.Scene;
import javafx.scene.control.Button;
import javafx.scene.layout.VBox;
import javafx.scene.web.WebView;
import javafx.stage.Modality;
import javafx.stage.Stage;
import netscape.javascript.JSObject;

public class CalibratedCPPInputEditor extends TreeDistributionInputEditor {

    List<String> taxa = Arrays.asList("Apple", "Banana", "Cherry");
    String chart = "";
    String d3 = "";

    private Stage popupStage;
    private WebView webView;
    
    private Object calibrations = null, idCounter = null;

    public CalibratedCPPInputEditor(BeautiDoc doc) {
        super(doc);
    }

    public CalibratedCPPInputEditor() {
		super();
	}

	@Override
    public Class<?> type() {
        return CalibratedCoalescentPointProcess.class;
    }

	@Override
	public void init(Input<?> input, BEASTInterface beastObject, int listItemNr, ExpandOption isExpandOption,
			boolean addButtons) {
		super.init(input, beastObject, listItemNr, isExpandOption, addButtons);

		try {
			loadjs();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
        CalibratedCoalescentPointProcess cpp = (CalibratedCoalescentPointProcess) m_beastObject;
        taxa = cpp.treeInput.get().getTaxonset().asStringList();

        Button calibrationButton = new Button("Manage calibrations");
        pane.getChildren().add(calibrationButton);
        calibrationButton.setOnAction(e -> openPopup());
	}
	
    private void loadjs() throws IOException {
    	{
	        InputStream instr = this.getClass().getClassLoader().getResourceAsStream("resources/chart.js");
	
	        // reading the files with buffered reader 
	        InputStreamReader strrd = new InputStreamReader(instr);
	       
	        BufferedReader rr = new BufferedReader(strrd);
	        String line;
	        // outputting each line of the file.
	        while ((line = rr.readLine()) != null) { 
	            //System.out.println(line);
	        	chart += line + "\n";
	        }
	        System.out.println("Chart " + chart.length());
    	}
    	{
	        InputStream instr = this.getClass().getClassLoader().getResourceAsStream("resources/d3.v7.min.js");
	    	
	        // reading the files with buffered reader 
	        InputStreamReader strrd = new InputStreamReader(instr);
	       
	        BufferedReader rr = new BufferedReader(strrd);
	        String line;
	        // outputting each line of the file.
	        while ((line = rr.readLine()) != null) { 
	            //System.out.println(line);
	        	d3 += line + "\n";
	        }
	        System.out.println("d3 " + d3.length());
    	}
    }
    
    private void openPopup() {
        popupStage = new Stage();
        popupStage.initModality(Modality.APPLICATION_MODAL);
        popupStage.setTitle("Edit List");

        webView = new WebView();

        
        // Bridge between JavaScript and Java
        webView.getEngine().getLoadWorker().stateProperty().addListener((obs, oldState, newState) -> {
            if (newState == javafx.concurrent.Worker.State.SUCCEEDED) {
                JSObject window = (JSObject) webView.getEngine().executeScript("window");
                window.setMember("javaConnector", new JavaConnector());
            }
        });
        
        String html = generateHtmlForm();
        webView.getEngine().loadContent(html);
        
        
        if (calibrations != null) {
        	//webView.getEngine().executeScript("updateNewick()");
        	webView.getEngine().getLoadWorker().stateProperty().addListener((obs, oldState, newState) -> {
                 if (newState == javafx.concurrent.Worker.State.SUCCEEDED) {
                 	webView.getEngine().executeScript("calibrations = JSON.parse('" + calibrations + "')");
                	webView.getEngine().executeScript("taxa = JSON.parse('"+ taxa +"')");
                	webView.getEngine().executeScript("idCounter = JSON.parse('"+ idCounter +"')");
                	
                	Object o =  webView.getEngine().executeScript("JSON.stringify(calibrations)");
                	System.out.println(o);
                	System.out.println(webView.getEngine().executeScript("JSON.stringify(taxa)"));
                	System.out.println(webView.getEngine().executeScript("JSON.stringify(idCounter)"));
                    // Call the function and get the return value
                     Object result = webView.getEngine().executeScript("renderCalibrations(); updateChart(); updateNewick();");
                     System.out.println("JavaScript function returned: " + result);
                 }
             });
       }
        
 
//        startPolling(webView);

        VBox vbox = FXUtils.newVBox();
        vbox.getChildren().add(webView);
        Button saveButton = new Button("Save");
        saveButton.setOnAction(e->getNewick());
        vbox.getChildren().add(saveButton);
        Scene popupScene = new Scene(vbox, 1024, 768);
        popupStage.setScene(popupScene);
        popupStage.showAndWait();
    }

	private String generateHtmlForm() {
        StringBuilder htmlBuilder = new StringBuilder();
        for (String item : taxa) {
            htmlBuilder.append(String.format(
                "\"%s\",", item
            ));
        }
        htmlBuilder.deleteCharAt(htmlBuilder.length()-1);
        String html = this.html;
        html = html.replaceFirst("let taxa = ", "let taxa = [" +htmlBuilder.toString() + "]");
        html = html.replaceFirst("</textarea>", "(" + htmlBuilder.toString().replaceAll("\\\"","")+ ")</textarea>");
        
        int i = html.indexOf("<SCRIPTS>");
        html = html.substring(0, i-1) +
        		"<script>" + chart + "</script>\n"
				+ "<script>" + d3 + "</script>\n"
				+ html.substring(i + 9); 

        return html;

//        StringBuilder htmlBuilder = new StringBuilder();
//        htmlBuilder.append("<html><body><h2>Edit List</h2><form id='myForm'>");
//        for (String item : items) {
//            htmlBuilder.append(String.format(
//                "<input type='text' name='item' value='%s' size='20'><br>", item
//            ));
//        }
//        htmlBuilder.append(
//            "<button type='button' onclick='submitForm()'>Save Changes</button>" +
//            "<script>" +
//            "function submitForm() {" +
//            "    var inputs = document.getElementsByName('item');" +
//            "    var result = [];" +
//            "    for (var i = 0; i < inputs.length; i++) {" +
//            "        result.push(inputs[i].value);" +
//            "    }" +
//            "    javaConnector.receiveList(JSON.stringify(result));" +
//            "    window.close();" +
//            "}" +
//            "</script>" +
//            "</form></body></html>"
//        );
//        return htmlBuilder.toString();
    }

    
    final String html = "<!DOCTYPE html>\n"
    		+ "<html lang=\"en\">\n"
    		+ "<head>\n"
    		+ "    <meta charset=\"UTF-8\">\n"
    		+ "    <meta name=\"viewport\" content=\"width=device-width, initial-scale=1.0\">\n"
    		+ "    <title>Calibration Lab</title>\n"
    		+ "<SCRIPTS>"
//    		+ "    <script src=\"https://cdn.jsdelivr.net/npm/chart.js\"></script>\n"
//    		+ "    <script src=\"https://d3js.org/d3.v7.min.js\"></script>\n"
    		+ "    <style>\n"
    		+ "        :root { --bg: #f8f9fa; --card: #ffffff; --primary: #2c3e50; --accent: #27ae60; --border: #dee2e6; --blue: #3498db; }\n"
    		+ "        body { font-family: 'Segoe UI', system-ui, sans-serif; background: var(--bg); margin: 0; color: var(--primary); font-size: 13px; }\n"
    		+ "        .app-layout { display: grid; grid-template-columns: 1fr 360px; height: 100vh; width: 100vw; overflow: hidden; }\n"
    		+ "        .col { display: flex; flex-direction: column; border-right: 1px solid var(--border); background: var(--card); height: 100vh; }\n"
    		+ "        .col-header { padding: 12px 15px; border-bottom: 1px solid var(--border); background: #fff; flex-shrink: 0; }\n"
    		+ "        .scroll-area { flex-grow: 1; overflow-y: auto; padding: 15px; background: #fdfdfd; scrollbar-width: thin; scrollbar-color: #ccc transparent; }\n"
    		+ "        .scroll-area::-webkit-scrollbar { width: 8px; }\n"
    		+ "        .scroll-area::-webkit-scrollbar-thumb { background-color: #cbd5e0; border-radius: 10px; border: 2px solid #fdfdfd; }\n"
    		+ "        \n"
    		+ "        /* Taxa Controls */\n"
    		+ "        .taxa-item { display: flex; gap: 5px; margin-bottom: 8px; }\n"
    		+ "        .taxa-input { flex-grow: 1; padding: 4px; border: 1px solid var(--border); border-radius: 4px; font-size: 12px; }\n"
    		+ "        .taxa-actions { padding: 12px; border-top: 1px solid #eee; display: flex; flex-direction: column; gap: 8px; }\n"
    		+ "\n"
    		+ "        /* Calibration Cards */\n"
    		+ "        .dist-card { background: #fff; border: 1px solid var(--border); border-radius: 8px; padding: 12px; margin-bottom: 15px; border-top: 5px solid #ccc; box-shadow: 0 2px 4px rgba(0,0,0,0.02); }\n"
    		+ "        .card-row { display: grid; grid-template-columns: 1fr 1fr; gap: 8px; margin: 10px 0; }\n"
    		+ "        .taxa-select-box { margin-top: 10px; padding: 8px; background: #f8f9fa; border-radius: 4px; border: 1px solid #eee; max-height: 180px; overflow-y: auto; }\n"
    		+ "        .selection-header { font-size: 10px; font-weight: bold; color: #444; margin: 8px 0 4px 0; border-bottom: 1px solid #ddd; padding-bottom: 2px; text-transform: uppercase; }\n"
    		+ "        .item-checkbox { display: flex; align-items: center; gap: 6px; font-size: 11px; margin-bottom: 4px; cursor: pointer; padding: 2px; border-radius: 3px; }\n"
    		+ "        .item-checkbox:hover { background: #eceff1; }\n"
    		+ "        .item-checkbox.disabled { color: #bbb; font-style: italic; cursor: not-allowed; opacity: 0.7; }\n"
    		+ "        .assigned-tag { color: #e74c3c; font-size: 9px; margin-left: auto; font-weight: bold; }\n"
    		+ "        \n"
    		+ "        #chart-container { height: 250px; margin-bottom: 15px; background: #fff; padding: 10px; border-radius: 8px; border: 1px solid #eee; flex-shrink: 0; }\n"
    		+ "        #tree-container { flex-grow: 1; background: #fff; position: relative; }\n"
    		+ "        \n"
    		+ "        .btn { padding: 6px 12px; border-radius: 4px; border: none; cursor: pointer; font-weight: 600; font-size: 11px; transition: 0.2s; text-align: center; }\n"
    		+ "        .btn-add { background: var(--accent); color: white; }\n"
    		+ "        .btn-fasta { background: var(--blue); color: white; }\n"
    		+ "        .btn-rm { background: #ff7675; color: white; padding: 2px 6px; border-radius: 3px; }\n"
    		+ "        \n"
    		+ "        h3 { margin: 0; font-size: 14px; text-transform: uppercase; color: #555; }\n"
    		+ "        label { font-size: 10px; font-weight: bold; color: #999; text-transform: uppercase; display: block; margin-bottom: 2px; }\n"
    		+ "        input[type=\"number\"], select { width: 100%; padding: 4px; border: 1px solid #ddd; border-radius: 3px; font-size: 12px; }\n"
    		+ "        .label-edit { border: none; border-bottom: 1px solid #eee; font-weight: bold; font-size: 13px; width: 80%; outline: none; }\n"
    		+ "    </style>\n"
    		+ "</head>\n"
    		+ "<body>\n"
    		+ "\n"
    		+ "<div class=\"app-layout\">\n"
    		+ "    <!-- LEFT: Taxa Manager -->\n"
    		+ "    <!--\n"
    		+ "    <div class=\"col\">\n"
    		+ "        <div class=\"col-header\"><h3>Taxa</h3></div>\n"
    		+ "        <div class=\"scroll-area\" id=\"taxaContainer\"></div>\n"
    		+ "        <div class=\"taxa-actions\">\n"
    		+ "            <button class=\"btn btn-add\" onclick=\"addTaxon()\">+ New Taxon</button>\n"
    		+ "            <button class=\"btn btn-fasta\" onclick=\"document.getElementById('fastaInput').click()\">📂 Load FASTA</button>\n"
    		+ "            <input type=\"file\" id=\"fastaInput\" style=\"display:none\" accept=\".fasta,.fa,.fas,.txt\" onchange=\"handleFasta(this)\">\n"
    		+ "        </div>\n"
    		+ "    </div>\n"
    		+ "    -->\n"
    		+ "\n"
    		+ "    <!-- MIDDLE: Calibrations -->\n"
    		+ "    <div class=\"col\">\n"
    		+ "        <div class=\"col-header\" style=\"display:flex; justify-content:space-between; align-items:center;\">\n"
    		+ "            <h3>Calibrations</h3>\n"
    		+ "            <button class=\"btn btn-add\" style=\"width:auto;\" onclick=\"addCalibration()\">+ Add Calibration</button>\n"
    		+ "        </div>\n"
    		+ "        <div id=\"chart-container\"><canvas id=\"mainChart\"></canvas></div>\n"
    		+ "        <div class=\"scroll-area\" id=\"middleScrollArea\">\n"
    		+ "            <div id=\"calibrationContainer\"></div>\n"
    		+ "        </div>\n"
    		+ "    </div>\n"
    		+ "\n"
    		+ "    <!-- RIGHT: Tree -->\n"
    		+ "    <div class=\"col\" style=\"border-right:none;\">\n"
    		+ "        <div class=\"col-header\"><h3>Phylogeny</h3></div>\n"
    		+ "        <div class=\"scroll-area\" style=\"padding:0; display:flex; flex-direction:column;\">\n"
    		+ "            <div style=\"padding:15px; border-bottom: 1px solid #eee;\">\n"
    		+ "                <label>Newick String</label>\n"
    		+ "                <textarea id=\"newickInput\" style=\"width:100%; height:40px; font-family:monospace; font-size:10px; padding:5px;\" readonly></textarea>\n"
    		+ "            </div>\n"
    		+ "            <div id=\"tree-container\"></div>\n"
    		+ "        </div>\n"
    		+ "    </div>\n"
    		+ "</div>\n"
    		+ "\n"
    		+ "<script>\n"
    		+ "    /** --- STATE --- **/\n"
    		+ "    let taxa = ;\n"
    		+ "    let calibrations = [];\n"
    		+ "    let idCounter = 0;\n"
    		+ "    let mainChart;\n"
    		+ "    const colors = ['#3498db', '#e74c3c', '#2ecc71', '#f1c40f', '#9b59b6'];\n"
    		+ "\n"
    		+ "    /** --- MATH ENGINE --- **/\n"
    		+ "    const DistMath = {\n"
    		+ "        gamma(z) {\n"
    		+ "            const p = [0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313, -176.61502916214059, 12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7];\n"
    		+ "            if (z < 0.5) return Math.PI / (Math.sin(Math.PI * z) * this.gamma(1 - z));\n"
    		+ "            z -= 1; let x = p[0]; for (let i = 1; i < 9; i++) x += p[i] / (z + i);\n"
    		+ "            let t = z + 7.5; return Math.sqrt(2 * Math.PI) * Math.pow(t, (z + 0.5)) * Math.exp(-t) * x;\n"
    		+ "        },\n"
    		+ "        beta(a, b) { return (this.gamma(a) * this.gamma(b)) / this.gamma(a + b); }\n"
    		+ "    };\n"
    		+ "\n"
    		+ "    const Configs = {\n"
    		+ "        bounded: { name: \"Bounded\", fields: [{id:'lower', label:'Lower', val:1.0}, {id:'upper', label:'upper', val:2.0}], pdf: (x, p) => x <= p.lower || x >= p.upper ? 0 : 1/(p.upper-p.lower), range: (p) => [p.lower, p.upper] },\n"
//    		+ "        lognormal: { name: \"Log-Normal\", fields: [{id:'m', label:'Log μ', val:1.0}, {id:'s', label:'Log σ', val:0.15}], pdf: (x, p) => x <= 0 ? 0 : (1/(x*p.s*Math.sqrt(2*Math.PI))) * Math.exp(-Math.pow(Math.log(x)-p.m,2)/(2*p.s*p.s)), range: (p) => [0, Math.exp(p.m + 3.5*p.s)] },\n"
//    		+ "        normal: { name: \"Normal\", fields: [{id:'m', label:'Mean', val:0}, {id:'s', label:'Sigma', val:1}], pdf: (x, p) => (1/(p.s*Math.sqrt(2*Math.PI))) * Math.exp(-Math.pow(x-p.m,2)/(2*p.s*p.s)), range: (p) => [p.m - 4*p.s, p.m + 4*p.s] },\n"
//    		+ "        gamma: { name: \"Gamma\", fields: [{id:'k', label:'Shape', val:2}, {id:'t', label:'Scale', val:2}], pdf: (x, p) => x <= 0 ? 0 : (1/(Math.pow(p.t, p.k)*DistMath.gamma(p.k)))*Math.pow(x, p.k-1)*Math.exp(-x/p.t), range: (p) => [0, (p.k*p.t)+(6*Math.sqrt(p.k)*p.t)] },\n"
//    		+ "        exponential: { name: \"Exp\", fields: [{id:'l', label:'Lambda', val:1}], pdf: (x, p) => x < 0 ? 0 : p.l * Math.exp(-p.l * x), range: (p) => [0, 5/p.l] },\n"
//    		+ "        beta: { name: \"Beta\", fields: [{id:'a', label:'Alpha', val:2}, {id:'b', label:'Beta', val:5}], pdf: (x, p) => (x<0 || x>1) ? 0 : (Math.pow(Math.max(1e-9, Math.min(1-1e-9, x)), p.a-1)*Math.pow(Math.max(1e-9, Math.min(1-1e-9, 1-x)), p.b-1))/DistMath.beta(p.a, p.b), range: (p) => [0, 1] },\n"
//    		+ "        cauchy: { name: \"Cauchy\", fields: [{id:'x0', label:'Loc', val:0}, {id:'g', label:'Scale', val:1}], pdf: (x, p) => (1/(Math.PI*p.g))*(1/(1+Math.pow((x-p.x0)/p.g, 2))), range: (p) => [p.x0-10*p.g, p.x0+10*p.g] },\n"
//   		+ "        laplace: { name: \"Laplace\", fields: [{id:'m', label:'Loc', val:0}, {id:'b', label:'Scale', val:1}], pdf: (x, p) => (1/(2*p.b))*Math.exp(-Math.abs(x-p.m) / p.b), range: (p) => [p.m-8*p.b, p.m+8*p.b] }\n"
    		+ "    };\n"
    		+ "\n"
    		+ "    /** --- FASTA IMPORT --- **/\n"
    		+ "    function handleFasta(input) {\n"
    		+ "        const file = input.files[0];\n"
    		+ "        if (!file) return;\n"
    		+ "\n"
    		+ "        const reader = new FileReader();\n"
    		+ "        reader.onload = function(e) {\n"
    		+ "            const text = e.target.result;\n"
    		+ "            const lines = text.split(/\\r?\\n/);\n"
    		+ "            const foundTaxa = [];\n"
    		+ "\n"
    		+ "            lines.forEach(line => {\n"
    		+ "                if (line.startsWith('>')) {\n"
    		+ "                    // Extract ID after '>' and before the first space\n"
    		+ "                    let name = line.substring(1).trim().split(/\\s+/)[0];\n"
    		+ "                    // Replace special chars for Newick safety\n"
    		+ "                    name = name.replace(/[^a-zA-Z0-9_]/g, '_');\n"
    		+ "                    if (name && !foundTaxa.includes(name)) {\n"
    		+ "                        foundTaxa.push(name);\n"
    		+ "                    }\n"
    		+ "                }\n"
    		+ "            });\n"
    		+ "\n"
    		+ "            if (foundTaxa.length > 0) {\n"
    		+ "                if (confirm(`Replace current taxa with ${foundTaxa.length} taxa found in FASTA?`)) {\n"
    		+ "                    taxa = foundTaxa;\n"
    		+ "                    // Reset calibrations assigned taxa because the names changed\n"
    		+ "                    calibrations.forEach(c => c.childrenTaxa = []);\n"
    		+ "                    renderTaxa();\n"
    		+ "                    renderCalibrations();\n"
    		+ "                    updateNewick();\n"
    		+ "                }\n"
    		+ "            } else {\n"
    		+ "                alert(\"No valid FASTA headers (starting with '>') were found.\");\n"
    		+ "            }\n"
    		+ "            input.value = ''; // Reset input\n"
    		+ "        };\n"
    		+ "        reader.readAsText(file);\n"
    		+ "    }\n"
    		+ "\n"
    		+ "    /** --- TAXA MANAGEMENT --- **/\n"
    		+ "    function addTaxon() { taxa.push(`Taxon_${taxa.length+1}`); renderTaxa(); renderCalibrations(); updateNewick(); }\n"
    		+ "    function removeTaxon(i) {\n"
    		+ "        const n = taxa[i]; taxa.splice(i, 1);\n"
    		+ "        calibrations.forEach(c => c.childrenTaxa = c.childrenTaxa.filter(t => t !== n));\n"
    		+ "        renderTaxa(); renderCalibrations(); updateNewick();\n"
    		+ "    }\n"
    		+ "    function updateTaxonName(i, name) {\n"
    		+ "        const old = taxa[i], n = name.replace(/\\s+/g, '_'); taxa[i] = n;\n"
    		+ "        calibrations.forEach(c => { const idx = c.childrenTaxa.indexOf(old); if(idx > -1) c.childrenTaxa[idx] = n; });\n"
    		+ "        renderCalibrations(); updateNewick();\n"
    		+ "    }\n"
    		+ "    function renderTaxa() {\n"
    		+ "        document.getElementById('taxaContainer').innerHTML = taxa.map((t, i) => `<div class=\"taxa-item\"><input type=\"text\" class=\"taxa-input\" value=\"${t}\" oninput=\"updateTaxonName(${i}, this.value)\"><button class=\"btn btn-rm\" onclick=\"removeTaxon(${i})\">✕</button></div>`).join('');\n"
    		+ "    }\n"
    		+ "\n"
    		+ "    /** --- CALIBRATIONS --- **/\n"
    		+ "    function addCalibration() {\n"
    		+ "        const id = idCounter++;\n"
    		+ "        const type = 'bounded';\n"
    		+ "        const params = {}; Configs[type].fields.forEach(f => params[f.id] = f.val);\n"
    		+ "        calibrations.push({ id, type, label: `Calibration_${id+1}`, color: colors[id % colors.length], childrenTaxa: [], childrenCals: [], params });\n"
    		+ "        renderCalibrations(); updateChart(); updateNewick();\n"
    		+ "    }\n"
    		+ "\n"
    		+ "    function isAncestor(childId, potentialParentId) {\n"
    		+ "        const child = calibrations.find(c => c.id === childId);\n"
    		+ "        if (!child) return false;\n"
    		+ "        if (child.childrenCals.includes(potentialParentId)) return true;\n"
    		+ "        return child.childrenCals.some(cid => isAncestor(cid, potentialParentId));\n"
    		+ "    }\n"
    		+ "\n"
    		+ "    function updateDistType(id, newType) {\n"
    		+ "        const cal = calibrations.find(c => c.id === id);\n"
    		+ "        cal.type = newType;\n"
    		+ "        cal.params = {}; Configs[newType].fields.forEach(f => cal.params[f.id] = f.val);\n"
    		+ "        renderCalibrations(); updateChart();\n"
    		+ "    }\n"
    		+ "\n"
    		+ "    function toggleExclusive(calId, type, item) {\n"
    		+ "        const currentCal = calibrations.find(c => c.id === calId);\n"
    		+ "        const listName = type === 'taxa' ? 'childrenTaxa' : 'childrenCals';\n"
    		+ "        const isSelected = currentCal[listName].includes(item);\n"
    		+ "        if (!isSelected) {\n"
    		+ "            calibrations.forEach(c => c[listName] = c[listName].filter(x => x !== item));\n"
    		+ "            currentCal[listName].push(item);\n"
    		+ "        } else { currentCal[listName] = currentCal[listName].filter(x => x !== item); }\n"
    		+ "        renderCalibrations(); updateChart(); updateNewick();\n"
    		+ "    }\n"
    		+ "\n"
    		+ "    function renderCalibrations() {\n"
    		+ "        const middlePane = document.getElementById('middleScrollArea');\n"
    		+ "        const mainScrollPos = middlePane ? middlePane.scrollTop : 0;\n"
    		+ "        const internalScrolls = {};\n"
    		+ "        document.querySelectorAll('.dist-card').forEach(card => {\n"
    		+ "            const calId = card.getAttribute('data-id');\n"
    		+ "            const box = card.querySelector('.taxa-select-box');\n"
    		+ "            if (box) internalScrolls[calId] = box.scrollTop;\n"
    		+ "        });\n"
    		+ "\n"
    		+ "        document.getElementById('calibrationContainer').innerHTML = calibrations.map(c => {\n"
    		+ "            const conf = Configs[c.type];\n"
    		+ "            return `\n"
    		+ "            <div class=\"dist-card\" data-id=\"${c.id}\" style=\"border-top-color: ${c.color}\">\n"
    		+ "                <div style=\"display:flex; justify-content:space-between; align-items:center;\">\n"
    		+ "                    <input type=\"text\" class=\"label-edit\" value=\"${c.label}\" oninput=\"calibrations.find(x=>x.id==${c.id}).label=this.value.replace(/\\\\s/g,'_'); updateChart(); updateNewick()\">\n"
    		+ "                    <button class=\"btn btn-rm\" onclick=\"calibrations=calibrations.filter(x=>x.id!=${c.id}); renderCalibrations(); updateChart(); updateNewick()\">✕</button>\n"
    		+ "                </div>\n"
    		+ "                <div style=\"margin-top:8px;\">\n"
    		+ "                    <label>Distribution Type</label>\n"
    		+ "                    <select onchange=\"updateDistType(${c.id}, this.value)\">\n"
    		+ "                        ${Object.keys(Configs).map(k => `<option value=\"${k}\" ${c.type === k ? 'selected' : ''}>${Configs[k].name}</option>`).join('')}\n"
    		+ "                    </select>\n"
    		+ "                </div>\n"
    		+ "                <div class=\"card-row\">\n"
    		+ "                    ${conf.fields.map(f => `<div><label>${f.label}</label><input type=\"number\" step=\"0.1\" value=\"${c.params[f.id]}\" oninput=\"calibrations.find(x=>x.id==${c.id}).params.${f.id}=parseFloat(this.value); updateChart()\"></div>`).join('')}\n"
    		+ "                </div>\n"
    		+ "                <div class=\"taxa-select-box\">\n"
    		+ "                    <div class=\"selection-header\">Child Taxa</div>\n"
    		+ "                    ${taxa.map(t => {\n"
    		+ "                        const owner = calibrations.find(o => o.id !== c.id && o.childrenTaxa.includes(t));\n"
    		+ "                        return `<label class=\"item-checkbox ${owner ? 'disabled' : ''}\"><input type=\"checkbox\" ${c.childrenTaxa.includes(t)?'checked':''} onchange=\"toggleExclusive(${c.id}, 'taxa', '${t}')\">${t} ${owner ? `<span class=\"assigned-tag\">${owner.label}</span>`:''}</label>`;\n"
    		+ "                    }).join('')}\n"
    		+ "                    <div class=\"selection-header\">Child Calibrations</div>\n"
    		+ "                    ${calibrations.filter(o => o.id !== c.id).map(o => {\n"
    		+ "                        const owner = calibrations.find(p => p.id !== c.id && p.childrenCals.includes(o.id));\n"
    		+ "                        const circular = isAncestor(o.id, c.id);\n"
    		+ "                        const disabled = (owner || circular) && !c.childrenCals.includes(o.id);\n"
    		+ "                        return `<label class=\"item-checkbox ${disabled ? 'disabled' : ''}\"><input type=\"checkbox\" ${c.childrenCals.includes(o.id)?'checked':''} ${disabled?'disabled':''} onchange=\"toggleExclusive(${c.id}, 'cals', ${o.id})\">${o.label} ${owner?`<span class=\"assigned-tag\">${owner.label}</span>`:''} ${circular?`<span class=\"assigned-tag\">LOOP</span>`:''}</label>`;\n"
    		+ "                    }).join('')}\n"
    		+ "                </div>\n"
    		+ "            </div>`;\n"
    		+ "        }).join('');\n"
    		+ "\n"
    		+ "        if (middlePane) middlePane.scrollTop = mainScrollPos;\n"
    		+ "        document.querySelectorAll('.dist-card').forEach(card => {\n"
    		+ "            const calId = card.getAttribute('data-id');\n"
    		+ "            const box = card.querySelector('.taxa-select-box');\n"
    		+ "            if (box && internalScrolls[calId] !== undefined) box.scrollTop = internalScrolls[calId];\n"
    		+ "        });\n"
    		+ "    }\n"
    		+ "\n"
    		+ "    /** --- VISUALIZATION --- **/\n"
    		+ "    function updateChart() {\n"
    		+ "        const datasets = calibrations.map(c => {\n"
    		+ "            const conf = Configs[c.type];\n"
    		+ "            const [min, max] = conf.range(c.params); const data = [];\n"
    		+ "            for(let i=0; i<=100; i++) { const x = min + (i * (max-min)/100); data.push({x, y: conf.pdf(x, c.params)}); }\n"
    		+ "            return { label: c.label, data, borderColor: c.color, backgroundColor: c.color+'22', fill: true, pointRadius: 0, tension: 0.3, showLine: true };\n"
    		+ "        });\n"
    		+ "        if (!mainChart) mainChart = new Chart(document.getElementById('mainChart'), { type: 'scatter', data: { datasets }, options: { responsive: true, maintainAspectRatio: false, animation: false, scales: { x: {type:'linear', title:{display:true, text:'Time'}}, y: {beginAtZero:true} } } });\n"
    		+ "        else { mainChart.data.datasets = datasets; mainChart.update(); }\n"
    		+ "    }\n"
    		+ "\n"
    		+ "    function updateNewick() {\n"
    		+ "        const hasParent = new Set();\n"
    		+ "        calibrations.forEach(c => { c.childrenTaxa.forEach(t => hasParent.add(t)); c.childrenCals.forEach(cid => hasParent.add(`CAL_${cid}`)); });\n"
    		+ "        const build = (cal) => {\n"
    		+ "            const els = [...cal.childrenTaxa, ...cal.childrenCals.map(id => {const ch = calibrations.find(x=>x.id===id); return ch ? build(ch) : '';})].filter(x => x!=='');\n"
    		+ "            return els.length ? `(${els.join(',')})${cal.label}` : cal.label;\n"
    		+ "        };\n"
    		+ "        const roots = [...calibrations.filter(c => !hasParent.has(`CAL_${c.id}`)).map(c => build(c)), ...taxa.filter(t => !hasParent.has(t))];\n"
    		+ "        const nwk = roots.length > 1 ? `(${roots.join(',')})Root;` : (roots[0] ? roots[0] + ';' : '');\n"
    		+ "        document.getElementById('newickInput').value = nwk;\n"
      		+ "        drawTree(nwk);\n"
//          + "    var result = [];" 
//          + "        result.push(10);" 
//          + "        result.push(12);" 
//          + "    javaConnector.receiveList(JSON.stringify(result));"
    		+ "        javaConnector.receiveList();" 	
    		+ "    }\n"
    		+ "\n"
    		+ "    function drawTree(newick) {\n"
    		+ "        const container = document.getElementById('tree-container'); container.innerHTML = '';\n"
    		+ "        if(!newick) return;\n"
    		+ "        try {\n"
    		+ "            const parse = (s) => {\n"
    		+ "                var e = [], r = {}, parts = s.split(/\\s*(;|\\(|\\)|,|:)\\s*/);\n"
    		+ "                for (var t = 0; t < parts.length; t++) {\n"
    		+ "                    var n = parts[t];\n"
    		+ "                    if (n === \"(\") { var c = {}; r.children = [c]; e.push(r); r = c; }\n"
    		+ "                    else if (n === \",\") { var c = {}; e[e.length-1].children.push(c); r = c; }\n"
    		+ "                    else if (n === \")\") { r = e.pop(); }\n"
    		+ "                    else if (n !== \":\" && n !== \";\" && parts[t-1] !== \":\") { r.name = n; }\n"
    		+ "                } return r;\n"
    		+ "            };\n"
    		+ "            const data = parse(newick);\n"
    		+ "            const width = container.clientWidth, height = container.clientHeight;\n"
    		+ "            const svg = d3.select(\"#tree-container\").append(\"svg\").attr(\"width\", width).attr(\"height\", height).append(\"g\").attr(\"transform\", \"translate(40,40)\");\n"
    		+ "            const root = d3.hierarchy(data); d3.cluster().size([height - 80, width - 150])(root);\n"
    		+ "            svg.selectAll(\"path\").data(root.links()).enter().append(\"path\").attr(\"fill\", \"none\").attr(\"stroke\", \"#ddd\").attr(\"stroke-width\", 2).attr(\"d\", d3.linkHorizontal().x(d => d.y).y(d => d.x));\n"
    		+ "            const nodes = svg.selectAll(\"g\").data(root.descendants()).enter().append(\"g\").attr(\"transform\", d => `translate(${d.y},${d.x})`);\n"
    		+ "            nodes.append(\"circle\").attr(\"r\", 5).attr(\"fill\", d => d.children ? \"#eee\" : \"#fff\").attr(\"stroke\", \"steelblue\").attr(\"stroke-width\", 2);\n"
    		+ "            nodes.append(\"text\").attr(\"dy\", \"0.31em\").attr(\"x\", d => d.children ? -12 : 12).style(\"text-anchor\", d => d.children ? \"end\" : \"start\").text(d => d.data.name || \"\").style(\"font-weight\", d => d.children ? \"bold\" : \"normal\");\n"
    		+ "        } catch(e) {}\n"
    		+ "    }\n"
    		+ "    \n"
    		+ "    renderTaxa(); addCalibration();\n"
    		+ "</script>\n"
//    		+ "<script type='text/javascript' src='https://web.archive.org/web/20130320104451/https://getfirebug.com/releases/lite/1.2/firebug-lite-compressed.js'></script>\n"
    		+ "</body>\n"
    		+ "</html>";
    
    private Object getNewick() {
		Object o = webView.getEngine().executeScript("document.getElementById('newickInput').value");
		System.out.println(o);
		calibrations = webView.getEngine().executeScript("JSON.stringify(calibrations)");
		// String taxa = (String) webView.getEngine().executeScript("JSON.stringify(taxa)");
		idCounter = webView.getEngine().executeScript("JSON.stringify(idCounter)");
		System.out.println(calibrations);
		return o;
	}

    // Inner class to handle JavaScript callbacks
    public class JavaConnector {
        public void receiveList() {
			Object o = webView.getEngine().executeScript("document.getElementById('newickInput').value");
			System.out.println(o);
        }
    }
}

