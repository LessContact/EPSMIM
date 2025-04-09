#ifndef ALIGN_HPP
#define ALIGN_HPP

#include <cstdlib>
#include <memory>
#include <new>

#include <cstdlib>
#include <memory>
#include <new>

template <typename T, std::size_t Alignment>
class AlignedAllocator {
public:
    using value_type = T;
    using pointer = T*;
    using const_pointer = const T*;
    using reference = T&;
    using const_reference = const T&;
    using size_type = std::size_t;
    using difference_type = std::ptrdiff_t;

    // Rebind allocator to type U
    template <typename U>
    struct rebind { using other = AlignedAllocator<U, Alignment>; };

    AlignedAllocator() = default;
    template <typename U>
    explicit AlignedAllocator(const AlignedAllocator<U, Alignment>&) {}

    pointer allocate(size_type n) {
        if (n == 0) {
            return nullptr;
        }

        if (Alignment < alignof(T)) {
            throw std::bad_alloc();  // Alignment must be at least the natural alignment of T
        }

        // std::aligned_alloc requires the size to be a multiple of alignment
        size_type alloc_size = n * sizeof(T);
        size_type padding = Alignment - (alloc_size % Alignment);
        if (padding != Alignment) {
            alloc_size += padding;
        }

        void* ptr = std::aligned_alloc(Alignment, alloc_size);
        if (ptr == nullptr) {
            throw std::bad_alloc();
        }

        return static_cast<pointer>(ptr);
    }

    void deallocate(pointer p, size_type) {
        std::free(p);
    }
};

template <typename T, std::size_t Alignment1, typename U, std::size_t Alignment2>
bool operator==(const AlignedAllocator<T, Alignment1>&, const AlignedAllocator<U, Alignment2>&) {
    return Alignment1 == Alignment2;
}

template <typename T, std::size_t Alignment1, typename U, std::size_t Alignment2>
bool operator!=(const AlignedAllocator<T, Alignment1>& lhs, const AlignedAllocator<U, Alignment2>& rhs) {
    return !(lhs == rhs);
}


#endif