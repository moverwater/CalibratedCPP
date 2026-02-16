package calibratedcpp.lphy.prior;

import lphy.base.evolution.tree.TimeTreeNode;

public class Calibration {
    String[] taxa;
    Double age;
    TimeTreeNode node;

    public Calibration(String[] taxa, Double age) {
        this.taxa = taxa;
        this.age = age;
    }

    public Calibration(String[] taxa) {
        this.taxa = taxa;
    }

    public Calibration(TimeTreeNode node) {
        this.node = node;
    }

    public String[] getTaxa() {
        return taxa;
    }

    public void setTaxa(String[] taxa) {
        this.taxa = taxa;
    }

    public Double getAge() {
        return age;
    }

    public void setAge(Double age) {
        this.age = age;
    }

    public TimeTreeNode getNode() {
        return node;
    }

    public void setNode(TimeTreeNode node) {
        this.node = node;
    }

    public String[] getTaxa(TimeTreeNode node) {
        return node.getLeafNames();
    }

}
