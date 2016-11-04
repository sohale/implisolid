#pragma once

#if AS_WORKER
#include <emscripten/emscripten.h>
#endif


extern "C" {
#if AS_WORKER
    // called from the main on front-end
    // All functions must have the same declaration: void(char* data, int size);
    void worker_api__started(char* data, int size);
    void worker_api__verts(char* data, int size);

    void workerapi_make_geometry(char* data, int size);
#endif
}

/* Worker API */
#if AS_WORKER
void worker_api__started(char* data, int size) {
    // int CHUNK_SIZE = 16;
    for(int i=0; i<10; i++) {
        std::cout << "Worker" << std::endl;
        emscripten_worker_respond_provisionally(0, 0);
    }
    emscripten_worker_respond(0, 0);
}
static vector<REAL> static_datum = {3.14, 1.4, 4 };
void worker_api__verts(char* data, int size) {

    auto ptr = get_v_ptr();
    auto sz = get_v_size() * sizeof(REAL) * 3;

    std::cout << "worker_api::verts()" << std::endl;
    // emscripten_worker_respond(static_cast<char*>(ptr), static_cast<int>(sz));
    emscripten_worker_respond(static_cast<char*>(static_cast<void*>(static_datum.data())), static_cast<int>(sizeof(static_datum)));
}

void workerapi_make_geometry_part1(char* data, int size) {
    //var data = {mp5:json, mcsettings:json}
}

#endif
/* End of Worker API */
