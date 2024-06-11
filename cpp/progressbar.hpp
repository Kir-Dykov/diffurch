#pragma once

#include <string>
#include <fstream>
#include <iostream>

class ProgressBar {
    public:
    
    int volume = 0;
    int progress_bar_length;
    string progress_bar_chars;
    string progress_chars = "`0123456789@";
    size_t progress_chars_size = progress_chars.size();
    int progress = 0;
    
    ProgressBar(int volume_) {
        volume = volume_;
        progress_bar_length = min(volume, 100);
        progress_bar_chars = string(progress_bar_length, '`');
    }
    
    
    void update() {
        double p = double(progress)/double(volume);
        double pp = p*progress_bar_length;
        if (int(pp) > 0) progress_bar_chars[int(pp)-1] = '@';
        if (int(pp) < progress_bar_length)
            progress_bar_chars[int(pp)]=progress_chars[int(progress_chars_size*(pp - int(pp)))];
    }
    
    
    void print(){
        cout << "\r" + progress_bar_chars << flush;
    }
};




