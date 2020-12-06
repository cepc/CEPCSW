#ifndef CLUPATRA_RUNTIMEMAP_H
#define CLUPATRA_RUNTIMEMAP_H
#include <map>

template<class U, class V>
class RuntimeMap {
    std::map<U, V> data;
    public:
    V& operator()(const U& u) {
        return data[u];
    }
    void clear() {
        data.clear();
    }
};
#endif
