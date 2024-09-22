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
    
    std::chrono::steady_clock::time_point begin;
    
    ProgressBar(int volume_) {
        volume = volume_;
        progress_bar_length = min(volume, 100);
        progress_bar_chars = string(progress_bar_length, '`');
        begin = std::chrono::steady_clock::now();
    }
    
    void increment() {
        progress++;
        update();
    }
    
    void update() {
        double p = double(progress)/double(volume);
        double pp = p*progress_bar_length;
        if (int(pp) > 0) progress_bar_chars[int(pp)-1] = '@';
        if (int(pp) < progress_bar_length)
            progress_bar_chars[int(pp)]=progress_chars[int(progress_chars_size*(pp - int(pp)))];
    }
    
    void print(){
        chrono::steady_clock::time_point end = chrono::steady_clock::now();
	    int seconds_passed = chrono::duration_cast<chrono::seconds>(end - begin).count();
        int seconds_total = int(seconds_passed * volume / progress);
        int seconds_left = seconds_total - seconds_passed;
        
        
        cout << "\r" + progress_bar_chars << " " 
            << (seconds_passed / 3600) << ":" << (seconds_passed / 60) % 60 << ":" << seconds_passed % 60 
            << " / "
            << (seconds_total  / 3600) << ":" << (seconds_total  / 60) % 60 << ":" << seconds_total  % 60 
            << "                                         " << flush;
    }
};




