#include "testApp.h"
//#include <armadillo>
//--------------------------------------------------------------
void testApp::setup(){

    ofSetFrameRate(500);

    float t_max;
    t_max = 0.1;

    n_osc = 200;
    fs = 5e4;

    n_t = t_max*fs;

    stimulus = new double[n_t+1];


    f0 = 1000;
    onset_dur = 0.005;

    for(int n = 0; n < n_t; n++)
    {
        t = n /fs;
        sig = 1*sin(6.28 * t * f0);

        if(t < onset_dur)
            sig =  sig*(1.0 - cos(PI * t / onset_dur)) / 2.0;

            stimulus[n] = sig;
    }

//    mat aux(stimulus,n_t+1,1,false);
  //  aux.save("stimulus.txt",raw_ascii);

    Y_t = new double[n_t*n_osc];
    V_t = new double[n_t*n_osc];

    n = 0;

    cochlea.setup(n_osc, n_t, fs, stimulus, Y_t, V_t);

    //Pass extra input parameters
    //cochlea.rho = rho;
    //cochlea.gaussian_elimination_init();


    //mat aux2(Y_t, n_osc,n_t, false);
    //aux2.save("Y_t.txt", raw_ascii);

    L = 600;
    dL = L/(float)n_osc;
    A = 1e0;
    zscale = 1;


}

//--------------------------------------------------------------
void testApp::update(){

    //Iterate
    if(n<n_t-10)
    {
        int n2 = n + 2;
        t = n2 /fs;
        float f0_t = ofMap(mouseX,0,ofGetWidth(),50,16000);
        f0 = f0 + (f0_t-f0)*0.001;

        sig = pow(6.28*f0,2.0)*sin(6.28 * t * f0);

        if(t < onset_dur)
            sig =  sig*(1.0 - cos(PI * t / onset_dur)) / 2.0;

            stimulus[n2] = sig;

        cochlea.update();
        n++;
    }

   // else
        //out << "fin" << endl;
    ofSetWindowTitle(ofToString(ofGetFrameRate()));
}

//--------------------------------------------------------------
void testApp::draw(){


    ofBackgroundGradient(ofColor(50,50,50), ofColor(0,0,0));

       ofSetColor(ofColor(255,255,255));

        ofNoFill();

        // draw the left channel:
        ofPushStyle();
            ofPushMatrix();
            ofTranslate(ofGetWidth()/2.0,ofGetHeight()/2.0,0);

            ofSetLineWidth(1);

                ofBeginShape();

                for (int i = 0; i < n_osc; i++)
                {

                    ofVertex(-L/2.0+(i+1)*dL,cochlea.Y[i]*A);
                }

                ofEndShape(false);

            ofPopMatrix();
        ofPopStyle();


}

//--------------------------------------------------------------
void testApp::keyPressed(int key){

}

//--------------------------------------------------------------
void testApp::keyReleased(int key){

}

//--------------------------------------------------------------
void testApp::mouseMoved(int x, int y){

}

//--------------------------------------------------------------
void testApp::mouseDragged(int x, int y, int button){

}

//--------------------------------------------------------------
void testApp::mousePressed(int x, int y, int button){

}

//--------------------------------------------------------------
void testApp::mouseReleased(int x, int y, int button){

}

//--------------------------------------------------------------
void testApp::windowResized(int w, int h){

}

//--------------------------------------------------------------
void testApp::gotMessage(ofMessage msg){

}

//--------------------------------------------------------------
void testApp::dragEvent(ofDragInfo dragInfo){

}
