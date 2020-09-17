#pragma once
#include <cassert>
#include <vector>

namespace Auction
{
namespace detail
{

/**
 * @brief class for easy locking of rows and cols
 *
 */
class Locks
{
  public:
    /**
     * @brief Create lock with defined size
     *
     * @param size number of rows/cols
     */
    Locks(size_t const & size)
        : _locks(size, false)
    {
    }

    Locks() = delete;

    /**
     * @brief locks the i-th flag
     *
     * @param index index of the flag
     */
    void lock(size_t const & index)
    {
        assert(index < _locks.size());
        if (!is_locked(index))
        {
            _locks[index] = true;
            _locked++;
        }
    }

    /**
     * @brief unlocks the i-th flag
     *
     * @param index index of the flag
     */
    void unlock(size_t const & index)
    {
        assert(index < _locks.size());

        if (is_locked(index))
        {
            _locks[index] = false;
            _locked--;
        }
    }

    /**
     * @brief indicates that the i-th flag is set locked
     *
     * @param index index = i-th flag
     * @return true if flag is set locked
     * @return false otherwise
     */
    bool is_locked(size_t const & index)
    {
        assert(index < _locks.size());
        return _locks[index];
    }

    /**
     * @brief indicates that all flags are locked
     *
     * @return true if all flags are locked
     * @return false otherwise
     */
    bool all_locked() const { return _locked == _locks.size(); }

    /**
     * @brief returns the number of the lock flags
     *
     * @return size_t number of the flags
     */
    size_t size() const { return _locks.size(); }

    /**
     * @brief returns the number of locked flags
     * 
     * @return size_t number of locked flags
     */
    size_t locked() const { return _locked; }

  private:
    std::vector<bool> _locks; ///< storage of the lock flags
    size_t _locked{0};        ///< number of locked entries
};

} // namespace detail
} // namespace Auction