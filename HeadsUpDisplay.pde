class HeadsUpDisplay {
    public boolean showFrameRate = true;
    public PFont font;
    public float fontSize;
    public color defaultColor = color(0, 0, 0);

    public HeadsUpDisplay() {
    }

    public void draw() {
        camera();
        hint(DISABLE_DEPTH_TEST);
        push();
        noLights();
        textMode(MODEL);
        textSize(24);
        textAlign(LEFT, TOP);
        text("fps "+round(frameRate), 10, 10);
        pop();
        hint(ENABLE_DEPTH_TEST);
    }

    public void setFont(PFont font) {
        this.font = font;
    }

    public void setFontSize(float fontSize) {
        this.fontSize = fontSize;
    }

    public void setColor(color c) {
        this.defaultColor = c;
    }
}