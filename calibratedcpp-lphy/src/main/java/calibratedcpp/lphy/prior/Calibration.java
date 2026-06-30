package calibratedcpp.lphy.prior;

import lphy.base.evolution.tree.TimeTreeNode;

public class Calibration {
    String[] taxa;
    Double age;
    TimeTreeNode node;
    Double upper;
    Double lower;

    public Calibration(String[] taxa, Double age) {
        this.taxa = taxa;
        this.age = age;
    }

    public Calibration(String[] taxa) {
        this.taxa = taxa;
    }

    public Calibration(String[] taxa, Double upper, Double lower) {
        this.taxa = taxa;
        this.upper = upper;
        this.lower = lower;
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

    public Double getUpper() {
        return upper;
    }

    public void setUpper(Double upper) {
        this.upper = upper;
    }

    public Double getLower() {
        return lower;
    }

    public void setLower(Double lower) {
        this.lower = lower;
    }

    public String[] getTaxa(TimeTreeNode node) {
        return node.getLeafNames();
    }

}
