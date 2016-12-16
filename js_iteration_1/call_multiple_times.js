
/**
    function callMultipleTimes (implisolid_, update_callback, burst_count, interval_msec);

    Asynchronously does a series of calls to the function  update_callback().
*/

'use strict';

var callMultipleTimes = function() {
    // private:
    var last_active = 0;
    var intervals_counter = 0;

    function function_body(burst_count, interval_msec, callback) {
        /** Asynchronously does a series of calls to the function  update_callback().
        Previous name: update_mc_mul tiple() */


        var BURST_COUNT = burst_count; //2; // 50;
        var INTERVAL_MSEC = interval_msec; //8-1; //40  // 1   // in Milliseconds
        // Works well on too small values of INTERVAL_MSEC (Semaphore).
        // The main algorithm performs much faster when this interval is larger (>100). Becasue of `gc`.


        var already_busy_doing_it=0;  // semaphore. (Refractory). Fast x1.
        var burst_counter = 0;

        intervals_counter++;
        if(intervals_counter>1) //dont start if there is already one running
        {
            intervals_counter--;
            return;
        }
        last_active = setInterval (
            function()
            {
                if(!last_active)
                    console.error("!last_active");

                burst_counter++; // mc_icounter

                console.error(BURST_COUNT);
                if(burst_counter<BURST_COUNT)  //This code is amazingly similar to Bursting dyamics in neurons.
                {

                    if(already_busy_doing_it<1) //removable
                    {
                        already_busy_doing_it++;

                        callback();

                        if(already_busy_doing_it>1)
                            console.error(">>>already_busy_doing_it:"+already_busy_doing_it);
                        if(!last_active)
                            console.error("!last_active");
                        already_busy_doing_it--;
                    }else{
                        console.log("hit");
                    }
                }
                else
                {
                    // The interval has finished its job. Time to go.
                    if(!last_active)
                        console.error("!last_active");
                    clearInterval(last_active); // But then don't we need to wait until the last one is finished?
                    intervals_counter--; // Allow next one in future
                    burst_counter = 0;
                    last_active = 0;
                }
            },
            INTERVAL_MSEC
        );
    }
    return function_body;
}();
