#include <Locks.hxx>
#include <gtest/gtest.h>

using namespace Auction::detail;

TEST(Locks, locks_has_predefined_size) 
{
    Locks locks(123);
    EXPECT_EQ(123, locks.size());
}

TEST(Locks, locking_increases_the_number_of_locked_flags)
{
    Locks locks(123);
    locks.lock(10);
    EXPECT_EQ(1, locks.locked());
}

TEST(Locks, unlocking_decreases_the_number_of_locked_flags)
{
    Locks locks(123);
    locks.lock(10);
    ASSERT_EQ(1, locks.locked());
    locks.unlock(10);
    EXPECT_EQ(0, locks.locked());
}

TEST(Locks, unlocking_the_same_flag_does_not_change_the_number_of_locked_flags)
{
    Locks locks(123);
    locks.lock(10);
    locks.lock(11);
    ASSERT_EQ(2, locks.locked());
    locks.unlock(10);
    EXPECT_EQ(1, locks.locked());
    locks.unlock(10);
    EXPECT_EQ(1, locks.locked());    
}

TEST(Locks, locking_the_same_flag_does_not_change_the_number_of_locked_flags)
{
    Locks locks(123);
    locks.lock(10);
    ASSERT_EQ(1, locks.locked());
    locks.lock(10);
    EXPECT_EQ(1, locks.locked());  
}

TEST(Locks, lock_flag_with_index_is_successful)
{
    Locks locks(123);
    locks.lock(10);
    EXPECT_TRUE(locks.is_locked(10));
}

TEST(Locks, unlock_flag_with_index_is_successful)
{
    Locks locks(123);
    locks.lock(10);
    ASSERT_TRUE(locks.is_locked(10));
    locks.unlock(10);
    EXPECT_FALSE(locks.is_locked(10));
}

TEST(Locks, locking_all_leads_to_all_locked_being_successful)
{
    Locks locks(5);
    for ( size_t i = 0; i < locks.size(); i++)
        locks.lock(i);
    EXPECT_TRUE(locks.all_locked());
}