/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package topicosespeciaisa;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author 11413771
 */
public class BioInformatica {

    private Set<String> arestas;
    private Set<String> vertices;

    public BioInformatica() {
        arestas = new HashSet<>();
        vertices = new HashSet<>();
    }

    public void calculaDistancia(String local) {
        try {
            FileReader arq = new FileReader(local);

            BufferedReader lerArq = new BufferedReader(arq);
            ArrayList<String> atomsE = new ArrayList<>();
            ArrayList<String> atomsI = new ArrayList<>();
            ArrayList<String> atomsLigantes = new ArrayList<>();
            ArrayList<String> hidrofobicos = new ArrayList<>();
            String linha = lerArq.readLine();
            Double distancia, atomEX, atomEY, atomEZ, atomIX, atomIY, atomIZ;

            hidrofobicos.add("ALA CB ");
            hidrofobicos.add("ARG CB");
            hidrofobicos.add("ARG CG");
            hidrofobicos.add("ARG CD");
            hidrofobicos.add("ARG CZ");
            hidrofobicos.add("ASN CB");
            hidrofobicos.add("ASN CG");
            hidrofobicos.add("ASP CB");
            hidrofobicos.add("ASP CG");
            hidrofobicos.add("CYS CB");
            hidrofobicos.add("GLN CB");
            hidrofobicos.add("GLN CG");
            hidrofobicos.add("GLN CD");
            hidrofobicos.add("GLU CB");
            hidrofobicos.add("GLU CG");
            hidrofobicos.add("GLU CD");
            hidrofobicos.add("HIS CB");
            hidrofobicos.add("HIS CG");
            hidrofobicos.add("HIS CD2");
            hidrofobicos.add("HIS CE1");
            hidrofobicos.add("ILE CB");
            hidrofobicos.add("ILE CG1");
            hidrofobicos.add("ILE CG2");
            hidrofobicos.add("ILE CD1");
            hidrofobicos.add("LEU CB");
            hidrofobicos.add("LEU CG");
            hidrofobicos.add("LEU CD1");
            hidrofobicos.add("LEU CD2");
            hidrofobicos.add("LYS CB");
            hidrofobicos.add("LYS CG");
            hidrofobicos.add("LYS CD");
            hidrofobicos.add("LYS CE");
            hidrofobicos.add("MET CB");
            hidrofobicos.add("MET CG");
            hidrofobicos.add("MET SD");
            hidrofobicos.add("MET CE");
            hidrofobicos.add("PHE CB");
            hidrofobicos.add("PHE CG");
            hidrofobicos.add("PHE CD1");
            hidrofobicos.add("PHE CD2");
            hidrofobicos.add("PHE CE1");
            hidrofobicos.add("PHE CE2");
            hidrofobicos.add("PHE CZ");
            hidrofobicos.add("PRO CB");
            hidrofobicos.add("PRO CG");
            hidrofobicos.add("PRO CD");
            hidrofobicos.add("SER CB");
            hidrofobicos.add("THR CB");
            hidrofobicos.add("THR CG2");
            hidrofobicos.add("TRP CG");
            hidrofobicos.add("TRP CD1");
            hidrofobicos.add("TRP CD2");
            hidrofobicos.add("TRP CE2");
            hidrofobicos.add("TRP CE3");
            hidrofobicos.add("TRP CZ2");
            hidrofobicos.add("TRP CZ3");
            hidrofobicos.add("TRP CH2");
            hidrofobicos.add("TYR CB");
            hidrofobicos.add("TYR CG");
            hidrofobicos.add("TYR CD1");
            hidrofobicos.add("TYR CD2");
            hidrofobicos.add("TYR CE1");
            hidrofobicos.add("TYR CE2");
            hidrofobicos.add("TYR CZ");
            hidrofobicos.add("VAL CB");
            hidrofobicos.add("VAL CG1");
            hidrofobicos.add("VAL CG2");

            while (linha != null) {

                if (linha.startsWith("ATOM") && linha.charAt(21) == 'E') {

                    atomsE.add(linha);

                }
                if (linha.startsWith("ATOM") && linha.charAt(21) == 'I') {

                    atomsI.add(linha);
                }
                linha = lerArq.readLine();
            }
            int i = 0;

            for (String atom : atomsE) {
                for (String hetatom : atomsI) {
                    if (hidrofobicos.contains(atom.substring(16, 20).trim() + " " + atom.substring(13, 16).trim())) {

                        atomEX = Double.parseDouble(atom.substring(30, 39));
                        atomEY = Double.parseDouble(atom.substring(38, 47));
                        atomEZ = Double.parseDouble(atom.substring(46, 55));
                        atomIX = Double.parseDouble(hetatom.substring(30, 39));
                        atomIY = Double.parseDouble(hetatom.substring(38, 47));
                        atomIZ = Double.parseDouble(hetatom.substring(46, 55));

                        distancia = Math.sqrt(Math.pow((atomEX - atomIX), 2) + Math.pow((atomEY - atomIY), 2) + Math.pow((atomEZ - atomIZ), 2));
                        if (distancia <= 7) {

                            String rotuloE = '"' + atom.substring(6, 11).trim() + "-" + atom.substring(16, 20).trim() + "-" + atom.substring(11, 15).trim() + '"';
                            String rotuloI = '"' + hetatom.substring(6, 11).trim() + "-" + hetatom.substring(16, 20).trim() + "-" + hetatom.substring(11, 15).trim() + '"';
                            String arestaEI = rotuloE + "|" + rotuloI;
                            vertices.add(rotuloE);
                            vertices.add(rotuloI);
                            arestas.add(arestaEI);
                        }
                    }
                }
            }
            System.out.println(vertices);
            System.out.println("---------");
            System.out.println(arestas);

            arq.close();
        } catch (IOException e) {
            System.err.printf("Erro na abertura do arquivo: %s.\n",
                    e.getMessage());
        }
    }

    public String parseEdges() {
        StringBuffer output = new StringBuffer();

        output.append("*Edges \n");
        List<String> verticesList = new ArrayList<>();
        verticesList.addAll(this.vertices);
        for (String aresta : this.arestas) {
            String nodes[] = aresta.split("\\|");
            int from = verticesList.indexOf(nodes[0]);
            int to = verticesList.indexOf(nodes[1]);
            output.append(from + "\t" + to + "\t1.000\n");

        }
        return output.toString();
    }

    public String parseVertices() {
        StringBuffer output = new StringBuffer();

        output.append("*Vertices " + this.vertices.size() + "\n");
        Iterator<String> iterator = this.vertices.iterator();
        int i = 1;
        while (iterator.hasNext()) {
            output.append(i + " " + iterator.next() + "\n");
            i++;
        }
        return output.toString();
    }

    public static void main(String[] args) throws FileNotFoundException {
        BioInformatica bioInformatica = new BioInformatica();
        PrintWriter writer = null;
        try {
           writer = new PrintWriter("tenten.net","UTF-8");
        } catch (UnsupportedEncodingException ex) {
            Logger.getLogger(BioInformatica.class.getName()).log(Level.SEVERE, null, ex);
        }
        String ppf = "1ppf.pdb";
        String r0r = "1r0r.pdb";
        bioInformatica.calculaDistancia(ppf);
        bioInformatica.parseEdges();
        writer.println(bioInformatica.parseVertices());
        writer.println("\n\n\n");
        writer.println(bioInformatica.parseEdges());
        writer.close();
        System.out.println();
    }

}
