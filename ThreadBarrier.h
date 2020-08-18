/*
 * ThreadBarrier.h
 *
 *  Created on: Jul 3, 2013
 *      Author: fb
 */

#ifndef THREADBARRIER_H_
#define THREADBARRIER_H_
#include <atomic>
#include <unistd.h>

namespace LSAP
{
/**
 * realizes a barrier with thread kill for std:threads
 */
class ThreadBarrier
{
public:
	/**
	 * creates a barrier which waits for n threads for continuation
	 * @param n number of threads
	 */
	ThreadBarrier(const size_t n, const bool & purespin = false, const size_t sleeptime = 5) :
			_n(n), _nwait(0), _step(0), _kill(false), _spin(purespin), _sleep(true), _sleeptime(sleeptime)
	{
	}
	virtual ~ThreadBarrier()
	{
	}

	void sleep() { if ( !_spin) _sleep = true; }

	void wakeup() { if ( !_spin) _sleep = false; }

	const bool operator()()
	{
		// return true if kill-signal
		if (_kill.load())
			return true;

		unsigned int step = _step.load();

		if (_nwait.fetch_add(1) == _n - 1)
		{
			_nwait.store(0);
			_step.fetch_add(1);
		}
		else
		{
			// waiting
			while (_step.load() == step)
			{
				if (_sleep) usleep(_sleeptime);
//				 return true if kill-signal
				if (_kill.load())
					return true;
			}
		}
		return false;
	}

	/**
	 * set kill flag for return value
	 */
	void kill()
	{
		_kill.store(true);
	}
protected:
	size_t _n;
	std::atomic<unsigned int> _nwait;
	std::atomic<unsigned int> _step;
	std::atomic<bool> _kill;
	bool _spin;
	bool _sleep;
	size_t _sleeptime;
};
} /* namespace LSAP */
#endif /* THREADBARRIER_H_ */
