class HeadsUpDisplay {
    public boolean showFrameRate;
    public PFont font;
    public float fontSize;
    public color defaultColor = color(0, 0, 0);

    public HeadsUpDisplay() {
    }

    public void draw() {
        push();
        noLights();
        PGraphics displayTexture = createGraphics(width,height);
        displayTexture.beginDraw();
        displayTexture.noStroke();
        displayTexture.textMode(MODEL);
        textAlign(LEFT,TOP);
        displayTexture.text("FPS: " + int(frameRate),50,70);
        displayTexture.endDraw();
        shader(unlitShader);
        beginShape(QUADS);
        texture(displayTexture);
        hint(DISABLE_DEPTH_TEST);
        vertex(0,0,0,0);
        vertex(width,0,displayTexture.width,0);
        vertex(width,height,displayTexture.width,displayTexture.height);
        vertex(0,height,0,displayTexture.height);
        endShape(CLOSE);
        pop();
        resetShader();
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